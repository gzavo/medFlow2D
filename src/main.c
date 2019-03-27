/*
 * Copyright (c) 2014, All Right Reserved, Gábor Závodszky, gabor@zavodszky.com

 * This source is subject to the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * Please see the License.txt file for more information.
 * All other rights reserved.

 * You should have received a copy of the license along with this
 * work in the License.txt file.
 * If not, see <http://creativecommons.org/licenses/by-nc-nd/4.0/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libgen.h>
#include "lbm/lb.h"
#include "lbm/boundaries.h"
#include "utils/arraylist.h"
#include "utils/queue.h"
#include "utils/utils.h"
#include "utils/slog.h"
#include "platelet/cell.h"
#include "inih/ini.h"

/*
#ifdef _OPENMP
	#include <omp.h>
#endif
*/

/*
 * Global variables
 */

Simulation sim;		// Global structure to hold simulation info

int LOG_LEVEL = 4;  // Above 2 its for debugging (set to 2 for normal operation)

// Main parameters of the simulation
real lbm_uMax;		// Maximum expected flow velocity in LBM units (should be <=0.1)
real phys_uMax;		// Maximum expected physical velocity based on input Re
real Re;			// Reynolds number
int L=0;			// Width of the inlet vessel in LBM units
real phys_nu;       // Physical viscosity
real lbm_nu;        // Lattice viscosity
real phys_rho;		// Fluid density in [kg/m^3]
real nu_ps_phys;
real nu_ps_lbm;
real phys_dx;			// Resolution -> size of one lattice in [m]
real rhoPSMean=1.0;		// Density for coupled scalar in LBM units

real defaultDensity = 1.0;		// Default LBM density -> 1.0 is good for numerical precision

int useMRT = 0;		// 0 to use bgk dynamics, 1 to use mrt

int useRegularized = 0; // 0 to use Zou-He, 1 to use regularised

int doPlatelets = 0;
int doAgeing = 0;

int numWorkers = 1;

// geometry description
char* geometryFileName;
PPMImage *geometry;

// boundary condition structures
char* inletScaleFileName;			// Text file that holds time dependent scale information
real outletDefaultDensity = 1.0;	// Default value, overwritten by XML data
ArrayList *inletNodes; 				// All nodes implementing velocityBC -> need to change inlet velocity during cardiac cycle
ArrayList *inletCoords;             // List of inlet numerical nodes
domainSide inletSide = UNSET;       // For automatic side selection
transientBC *inletBC;				// Transient data structure
int calcParabolicProfile =0;        // Whether to compute parabolic profile for inlet
real Uin_x = 0;                     // Inlet velocity directions (normalised!)
real Uin_y=0;
int autoDir = 0;                    // whether to calculate normal direction or used given one

// Warm up sequence values
real maxULBMchangePerSec = 0.1;		// Largest change during (will be divided for iterations)

// run options
real runTime;						// total run time of the simulation [s]
int monitorUpdateInterval = 25;		// Update current image every 25 iteration
real saveInterval = 0.01;		    // Save output freq. in [s]
char *dirBase;					    // Base output directory for saving images
char dirMonitor[1024];				// Base directory for saving current simulation state information
char dirResData[1024];              // Base directory for saving data files
char dirResImg[1024];               // Base directory to save image files
char *fileLog = "medFlow2D_solver";      // Log file
char *fileCheckpoint = "checkpoint.dat"; // File to sava to, or load from the checkpoint
int doSaveImg =0;					// Whether we want images to be saved
int doSaveData =0;					// Whether we want full saving or just images
int doSaveCoag=0;                   // Whether to save near wall data along with other save files
int useCompression=1;               // Whether to use ZIP compression for output files
int saveSparseStep=10;				// When saving sparse data step this big in every direction
int phase=0;						// Global simulation phase (0 - initialization, >0 simulating)
int useCheckPoint=0; 				// Whether to use checkpoint file
char *fileInputINI;                 // Input ini file
char *dirInput;                     // Input file directory
char *fullCmdFile;

// TODO: Account for this settings during initialisation and destruction
int isCoupled = 0;					// Whether to calculate with platelets or not
bool computeStressTensor = false;	// Compute the stress tensor
bool useSmagorinsky = false;		// Whether to use turbulence modeling
real smagorinskyC = 0.038;          // Smagorinsky constant
real margForceRatio = 1.0;


int glycocalyxWidth = 0; 			// Number of cells to omit forcedBgk next to walls. It can account for glycocalyx layers while also improving the stability

// non-Newtonian dynamics
int useNonNewt = 0;					// Whether to use non-Newtonian dynamics for fluid material modeling
real CYlambda=0;
real CYa=0;
real CYn=0;
real CYnu_0=0;
real CYnu_inf=0;

// Thrombus modelling
int calcThrombus = 0;				// Whether to use thrombus model (couple simulation required)
ArrayList *nearWallCells;			// Collection of cells that can turn into coagwall
int coagulationUpdateInterval = 0;	// Update near wall cell history requency
real activationTime = 2e-2;			// Longest required activation time [s] (like ADP = 20ms)
ArrayList *woundCoords;				// Store starting wound positions
real woundMultiplier=1;				// Increase coag probability at wounded cells

// Porous material options
int calcPorousConsts =0;			// Caclulate porous constants based on sphere packed bed analityc solution
real c0, c1;						// Porous constant for Darcy term and Forcheimer term respectively
real porosity, diameter;			// Values for porous calculations

// output image color scaling
real imgRhoScaleMin = 0.85;
real imgRhoScaleMax = 1.15;
real imgPSScaleMin = 0.8;
real imgPSScaleMax = 1.2;
real imgRhoScale = 0.3;
real imgPSScale = 0.4;

// Generate full file path
char *filePath(const char * file, const char* dir )
{
	char *f = (char*)calloc(256, sizeof(char));
	strcpy(f, dir); strcat(f, "/");
	strcat(f, file);
	return f;
}

// Calculate analytical porous constants based on spherical particle bed
void calcPorousConstants(Simulation *sim, real nuPhys, real porosity, real diameter, real *c0, real *c1)
{
	real Fe, K, eps, dp;

	eps = porosity;//0.8;				// porosity
	dp = diameter;//0.001;				// particle diameter (controls geometric function)
	K = eps*eps*eps*dp*dp/(150.0*(1.0-eps)*(1.0-eps));	// permeability
	Fe = 1.75/(sqrt(150.0*eps*eps*eps));		// geometric function

		/*
		c0 = 0.5*(1+eps*0.5*(nuPhys/K));		// Linear part of force
		c1 = eps*0.5*Fe/sqrt(K);			// Quadratic part of force
		*/
	*c0 = eps*nuPhys/K * sim->dt;			// Linear part of force
	*c1 = eps*Fe/sqrt(K) * sim->dx;			// Quadratic part of force

	//F[0] = omega; F[1]=c0; F[2]=c1;
}

// Process configuration file entries
int handler(void* user, const char* section, const char* name, const char* value) {
	//configuration* pconfig = (configuration*)user;

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
	// PARAMETERS
	if (MATCH("parameters", "uMax_lbm")) {
		lbm_uMax = atof(value);
	} else if (MATCH("parameters", "viscosity")) {
		phys_nu = atof(value);
	} else if (MATCH("parameters", "density")) {
		phys_rho = atof(value);
	} else if (MATCH("parameters", "xRes")) {
		phys_dx = atof(value);
	} else if (MATCH("parameters", "useCheckPoint")) {
		useCheckPoint = atoi(value);
	} else if (MATCH("parameters", "checkPointFile")) {
		fileCheckpoint = strdup(value);
	} else if (MATCH("parameters", "useNonNewtonian")) {
		useNonNewt = atoi(value);
	} else if (MATCH("parameters", "useSmagorinsky")) {
		useSmagorinsky = (bool) atoi(value);
	} else if (MATCH("parameters", "smagorinskyC")) {
        smagorinskyC = (real) atof(value);
    }  else if (MATCH("parameters", "useMRT")) {
		useMRT = atoi(value);
	} else if (MATCH("parameters", "useRegularized")) {
		useRegularized = atoi(value);
	} // COUPLED FIELD
	else if (MATCH("coupled_field", "rhoMean")) {
		rhoPSMean = atof(value);
	} else if (MATCH("coupled_field", "psDiffusion")) {
		nu_ps_phys = atof(value);
	} else if (MATCH("coupled_field", "coupled")) {
		isCoupled = atoi(value);
	} else if (MATCH("coupled_field", "doAgeing")) {
        doAgeing = atoi(value);
    } else if (MATCH("coupled_field", "doPlatelets")) {
        doPlatelets = atoi(value);
    } else if (MATCH("coupled_field", "margForceRatio")) {
		margForceRatio = atof(value);
	} else if (MATCH("coupled_field", "glycocalyxWidth")) {
		glycocalyxWidth = atoi(value);
	} // THROMBUS
    else if (MATCH("thrombus", "calcThrombus")) {
        calcThrombus = atoi(value);
    } else if (MATCH("thrombus", "activationTime")) {
    	activationTime = atof(value);
    } else if (MATCH("thrombus", "maxHistory")) {
    	setMaxHistory(atoi(value));
    } else if (MATCH("thrombus", "rhoADPDecr")) {
		setRhoADPDecrement(atof(value));
    } else if (MATCH("thrombus", "defRhoADP")) {
		setDefaultRhoADP(atof(value));
	} else if(MATCH("thrombus", "coagCoeff")) {
		setCoagCoeff(atof(value));
    } else if (MATCH("thrombus", "woundMultiplier")) {
    	woundMultiplier = atof(value);
    } // GEOMETRY
	else if (MATCH("geometry", "file")) {
        geometryFileName = strdup(value);
		geometry = readPPM(filePath(geometryFileName, dirInput));
    } // RUN OPTIONS
	else if (MATCH("run_config", "run_time")) {
        runTime = atof(value);
    } else if (MATCH("run_config", "workers")) {
		numWorkers = atoi(value);
	} // BOUNDARY OPTIONS
	else if (MATCH("boundary", "inlet_scale")) {
        inletScaleFileName = strdup(value);
    } else if (MATCH("boundary", "outlet_density")) {
        outletDefaultDensity = atof(value);
    } else if (MATCH("boundary", "Re")) {
	    Re = atof(value);
    }  else if (MATCH("boundary", "parabolic_profile")) {
		calcParabolicProfile = atoi(value);
	} else if (MATCH("boundary", "uin_x")) {
		Uin_x = atof(value);
	} else if (MATCH("boundary", "uin_y")) {
		Uin_y = atof(value);
	} // POROUS
    else if (MATCH("porous", "calcPorousConsts")) {
    	calcPorousConsts =atoi(value);
    } else if (MATCH("porous", "porosity")) {
    	porosity=atof(value);
    } else if (MATCH("porous", "diameter")) {
    	diameter=atof(value);
    } else if (MATCH("porous", "c0")) {
    	c0=atof(value);
    } else if (MATCH("porous", "c1")) {
    	c1=atof(value);
    } // WARMUP
    else if (MATCH("warmup", "maxULBMchangePerSec")) {
    	maxULBMchangePerSec=atof(value);
    } // OUTPUT
    else if(MATCH("output", "saveInterval")) {
    	saveInterval = atof(value);
    } else if(MATCH("output", "workingDir")) {
    	dirBase = strdup(value);
    } else if(MATCH("output", "save_images")) {
        doSaveImg = atoi(value);
    } else if(MATCH("output", "save_coagdata")) {
        doSaveCoag = atoi(value);
    } else if(MATCH("output", "save_datafiles")) {
    	doSaveData = atoi(value);
    }else if(MATCH("output", "useCompression")) {
		useCompression = atoi(value);
	} else if (MATCH("output", "monitorUpdateInterval")) {
		monitorUpdateInterval =atoi(value);
	} else if(MATCH("output", "imgRhoScale")) {
    	imgRhoScale = atof(value);
    } else if(MATCH("output", "imgPSScale")) {
    	imgPSScale = atof(value);
    } else if(MATCH("output", "saveSparseStep")) {
    	saveSparseStep = atoi(value);
    } // Carreau-Yasuda
    else if(MATCH("carreau-yasuda", "nu_0")) {
        CYnu_0 = atof(value);
    } else if(MATCH("carreau-yasuda", "nu_inf")) {
    	CYnu_inf = atof(value);
    } else if(MATCH("carreau-yasuda", "a")) {
    	CYa = atof(value);
    } else if(MATCH("carreau-yasuda", "lambda")) {
    	CYlambda = atof(value);
    } else if(MATCH("carreau-yasuda", "n")) {
    	CYn = atof(value);
    }

	// TODO - More parameters ?! 	
	//else { return 0;  /* unknown section/name, error */    }
	
    return 1;
}

// Wrapper for configuration file parsing
int loadSetupIniFile(char *fileName) {

	slog(1, SLOG_INFO, "-> Loading simulation ini file...\n");
		
	int error;
	
	error = ini_parse(fileName, handler, NULL);
    if (error < 0) {
        slog(0, SLOG_ERROR, "Can't read '%s'!\n", fileName);
        return 2;
    }
    else if (error) {
        slog(0, SLOG_ERROR, "Bad configuration file (first error on line %d)!\n", error);
        return 3;
    }
	
	slog(1, SLOG_INFO, "-> Loading done successfully!\n");
	
	return 0;
}

// Use this during simulation to set boundary conditions
// Warning: this assumes small changes, and does not change the underlying distribution to equilibrium!
void setVelocityBoundary(ArrayList *dyn, int length, domainSide side, real uMax_lbm)
{
	int i;
	for (i = 0; i < length; i++) {
		Dynamics *d = (Dynamics *) dyn->data[i];
		VelocityBCData *bc = (VelocityBCData *) d->selfData;

		real vel;

		// Select profile type
		if(calcParabolicProfile) {
			int R = floor(length / 2);
			int r = i - R;
			vel = poiseuilleProfile(r, uMax_lbm, R);
		}
		else
			vel = uMax_lbm;

		// Set velocity direction
		if(autoDir) {
			switch (side) {
				case TOP:
					bc->ux = 0;
			        bc->uy = vel;
			        break;
				case BOTTOM:
					bc->ux = 0;
			        bc->uy = -vel;
			        break;
				case LEFT:
					bc->ux = vel;
			        bc->uy = 0;
			        break;
				case RIGHT:
					bc->ux = -vel;
			        bc->uy = 0;
			        break;
				case UNSET:
					slog(1, SLOG_ERROR, "Error setting velocity boundary: Unset side!!!\n");
			        break;
			}
		}
		else
		{
			bc->ux = Uin_x * vel;
			bc->uy = Uin_y * vel;
		}
	}
}

// Allocate dynamics and values to every cell independently -> parallelization
int constructSimulation()
{
	slog(1, SLOG_INFO, "-> Allocating memory for grid (%d, %d)\n", geometry->x, geometry->y);
	
	// For now, always allocate a coupled simulation (if you don't like it use Palabos...)
	// Stores [0..x+1] x [0..y+1] nodes, for later ghost nodes
	constructSim(&sim, geometry->x, geometry->y, 0, isCoupled, calcThrombus);

    if(numWorkers > 1) {
        slog(1, SLOG_INFO, "-> Parallel support is enabled!\n");
        slog(1, SLOG_INFO, "-> Number of worker threads: %d\n", numWorkers);
    }

/*	int hasOMP=0;

	#ifdef _OPENMP
    	//Set the number of worker threads
    	omp_set_num_threads(numWorkers);
		hasOMP=1;
	#endif

	slog(1, SLOG_INFO, "-> OpenMP support enabled: %d\n", hasOMP);
*/

	slog(1, SLOG_INFO, "-> Computing simulation parameters...\n");

	// Check inlet directions
	if(Uin_x == 0.0 && Uin_y == 0.0)
	{
		slog(1, SLOG_INFO, "-> Using automatic normal inlet flow.\n");
		autoDir = 1;
	}
	else
	{
		real inMag = sqrt(Uin_x*Uin_x + Uin_y*Uin_y);
		Uin_x /= inMag;
		Uin_y /= inMag;
		slog(1, SLOG_INFO, "-> Using normalised inlet direction: (%f, %f).\n", Uin_x, Uin_y);
		autoDir = 0;
	}

	// Get inlet width (only one inlet is allowed in the current implementation!)
	int i,j;
	for(i=0; i<geometry->x; i++)
		for(j=0; j<geometry->y; j++) {
			int r,g,b;
			int idx = j*geometry->x+i;
			r = geometry->data[idx].r;
			g = geometry->data[idx].g;
			b = geometry->data[idx].b;

			if(r==255 && g==0 && b==0)
                L++; // Count the width of the inlet opening
		}

    slog(1, SLOG_INFO, "-> Inlet width: %d numeric cell.\n", L);

	sim.Re = Re;
	sim.dx = phys_dx;

    if(calcParabolicProfile) {
        slog(1, SLOG_INFO, "-> Using parabolic inlet profile...\n");
		lbm_nu = ((2.0 / 3.0) * lbm_uMax) * L / Re;    // kinematic fluid viscosity in LBM units, in case of a parabolic profile
	}
	else {
        slog(1, SLOG_INFO, "-> Using bottleneck inelet profile.\n");
        lbm_nu = (lbm_uMax) * L / Re;                // kinematic fluid viscosity in LBM units, in case of a block profile
    }

	// Parameter sanity check
    if(phys_dx == 0 || Re == 0 || L == 0 || lbm_uMax == 0)
	{
        slog(1, SLOG_ERROR, "Error during parameters sanity check!\nSomething might be zero that is not supposed to.\n");
        freeSimulation();
        exit(-1);
	}

	sim.omega = 1. / (3 * lbm_nu + 1./2.);      		// relaxation parameter
	sim.dt = phys_dx * phys_dx * lbm_nu / phys_nu;				// Get physical time step

	//sim.dm = phys_rho*sim.dx*sim.dx*sim.dx/defaultDensity;
    sim.dm = phys_rho * sim.dx * sim.dx;	// Get unit density in a cell

    phys_uMax = lbm_uMax * sim.dx / sim.dt;			// Maximum possible velocity

	real nu_ps_lbm = nu_ps_phys/ phys_dx / phys_dx *sim.dt;	// diffusion coefficient
	sim.coupledOmega = 1. / (3*nu_ps_lbm+0.5);		// omega for the coupled passive scalar field

	// TODO: Do we need this scaling at all?
	real defRhoADP = getDefaultRhoADP()*sim.dx*sim.dx; // [kg/lattice]
	//setDefaultRhoADP(defRhoADP);

	real coagCoeff = getCoagCoeff()/sim.dx/sim.dx;	// [tau_max/m^3]
	//setCoagCoeff(coagCoeff);

	//TODO: Not really here but create c0,c1 for the appropriate material width -> stent width
	// Typical stent strut is between 65-115e-6 m

	// Change c0, c1 from SI to LBM units
	// c0 = mu/K = [Pa*s]/[m^2] = [kg/m^3*s^2]

	//c0 = c0*sim.dx*sim.dx*sim.dx*sim.dt*sim.dt/sim.dm;
    c0 = c0*sim.dx*sim.dx*sim.dt*sim.dt/sim.dm;

	// c1 = [rho*c_f/sqrt(K)] = [kg/m^3]/[m] = [kg/m^4]

	//c1 = c1*sim.dx*sim.dx*sim.dx*sim.dx/sim.dm;
    c1 = c1*sim.dx*sim.dx*sim.dx/sim.dm;

	// non-Newtonian stuff
	if(useNonNewt)
	{
		sim.useNonNewtonian = 1;

		// Non-dimensional parameters
		sim.nonNewt.n = CYn;
		sim.nonNewt.a = CYa;

		// Dimensional parameters, change to LBM units
		sim.nonNewt.lambda 	= CYlambda/sim.dt;
		sim.nonNewt.nu_0	= CYnu_0*sim.dt/sim.dx/sim.dx;
		sim.nonNewt.nu_inf 	= CYnu_inf*sim.dt/sim.dx/sim.dx;
	}

    sim.useSmagorinsky = useSmagorinsky;

	if(doPlatelets || useNonNewt || useSmagorinsky) {
		// Enable shear stress computation - needed for drift force calculations
		computeStressTensor = true;	// this forces shear stress computation at every fluid node -> platelet macroscopic model needs it
	}

    imgRhoScaleMin = 1.0 - imgRhoScale/2.0;
    imgRhoScaleMax = 1.0 + imgRhoScale/2.0;

    if(doAgeing) {
        imgPSScaleMin = 0;
        imgPSScaleMax = imgPSScale;
    }
    else {
        imgPSScaleMin = rhoPSMean - imgPSScale/2.0;
        imgPSScaleMax = rhoPSMean + imgPSScale/2.0;
    }

	// TODO: do some sanity check here for values, warn if something is extremely out of order

	slog(1, SLOG_INFO, "-> Choosing dynamics...\n");

	// Main field dynamics

	void (*dynamicsFun) (real* fPop, void* selfData, NodeStore* store);

	void (*leftVelBC) (real* fPop, void* selfData, NodeStore* store);
	void (*rightVelBC) (real* fPop, void* selfData, NodeStore* store);
	void (*topVelBC) (real* fPop, void* selfData, NodeStore* store);
	void (*bottomVelBC) (real* fPop, void* selfData, NodeStore* store);

	void (*leftPresBC) (real* fPop, void* selfData, NodeStore* store);
	void (*rightPresBC) (real* fPop, void* selfData, NodeStore* store);
	void (*topPresBC) (real* fPop, void* selfData, NodeStore* store);
	void (*bottomPresBC) (real* fPop, void* selfData, NodeStore* store);


	if(useMRT) {
	    slog(1, SLOG_INFO, "\t-> Using MRT dynamics.\n");
	    dynamicsFun = &mrt;
    }
    else {
	    slog(1, SLOG_INFO, "\t-> Using BGK dynamics.\n");
	    dynamicsFun = &bgk;
    }

    if(useSmagorinsky)
        slog(1, SLOG_INFO, "\t-> Using Smagorinsky turbulence model (C=%f).\n", smagorinskyC);
    else
        slog(1, SLOG_INFO, "\t-> Using laminar model.\n");

	if(useRegularized) {
		slog(1, SLOG_INFO, "\t-> Using regularized BCs.\n");
		leftVelBC = &leftRegularized;
		rightVelBC = &rightRegularized;
		topVelBC = &topRegularized;
		bottomVelBC = &bottomRegularized;

		leftPresBC = &leftPressureRegularized;
		rightPresBC = &rightPressureRegularized;
		topPresBC = &topPressureRegularized;
		bottomPresBC = &bottomPressureRegularized;
	}
	else {
		slog(1, SLOG_INFO, "\t-> Using Zou-He BCs.\n");
		leftVelBC = &leftZouHe;
		rightVelBC = &rightZouHe;
		topVelBC = &topZouHe;
		bottomVelBC = &bottomZouHe;

		leftPresBC = &leftPressureZouHe;
		rightPresBC = &rightPressureZouHe;
		topPresBC = &topPressureZouHe;
		bottomPresBC = &bottomPressureZouHe;
	}

    // Coupled field dynamics
    void (*dynamicsFunPS) (real* fPop, void* selfData, NodeStore* store);
    void (*dynamicsFunPSProt) (real* fPop, void* selfData, NodeStore* store);

    void (*leftVelBCPS) (real* fPop, void* selfData, NodeStore* store);
    void (*rightVelBCPS) (real* fPop, void* selfData, NodeStore* store);
    void (*topVelBCPS) (real* fPop, void* selfData, NodeStore* store);
    void (*bottomVelBCPS) (real* fPop, void* selfData, NodeStore* store);

    void (*leftPresBCPS) (real* fPop, void* selfData, NodeStore* store);
    void (*rightPresBCPS) (real* fPop, void* selfData, NodeStore* store);
    void (*topPresBCPS) (real* fPop, void* selfData, NodeStore* store);
    void (*bottomPresBCPS) (real* fPop, void* selfData, NodeStore* store);


    // Check for thrombus model requirements
    if(calcThrombus)
        if(!isCoupled || !doPlatelets) {
            slog(0, SLOG_ERROR, "ERROR: the thrombus model requires a coupled simulation with platelet modeling!\n");
            freeSimulation();
            exit(-1);
        }

    // Check advanced dynamics requirements
    if(isCoupled)
        if(doAgeing && doPlatelets){
            slog(0, SLOG_ERROR, "ERROR: Only one coupled dynamics can be active at a time!\n");
            freeSimulation();
            exit(-1);
        }
        else if (!doAgeing && !doPlatelets) {
            slog(0, SLOG_ERROR, "ERROR: Requested a coupled simulation, however, none of the coupled dynamics is enabled!\n");
            freeSimulation();
            exit(-1);
        }

    if(doPlatelets) {
        slog(1, SLOG_INFO, "\t-> Applying platelet dynamics.\n");

        dynamicsFunPS = &bgkForcedPS;
        dynamicsFunPSProt = &bgkPS;

        // No velocity boundary for passive scalar
        leftVelBCPS = &leftPressurePS;
        rightVelBCPS = &rightPressurePS;
        topVelBCPS = &topPressurePS;
        bottomVelBCPS = &bottomPressurePS;

        leftPresBCPS = &leftPressurePS;
        rightPresBCPS = &rightPressurePS;
        topPresBCPS = &topPressurePS;
        bottomPresBCPS = &bottomPressurePS;
    }
    else {//if (doAgeing){
        slog(1, SLOG_INFO, "\t-> Applying fluid ageing dynamics.\n");

        rhoPSMean = 0.0; // Start ageing from 0

        dynamicsFunPS = &ageingDyn;
        dynamicsFunPSProt = &ageingDyn;

        leftVelBCPS = &zeroAge;
        rightVelBCPS = &zeroAge;
        topVelBCPS = &zeroAge;
        bottomVelBCPS = &zeroAge;

        leftPresBCPS = &leftAgeOutlet;
        rightPresBCPS = &rightAgeOutlet;
        topPresBCPS = &topAgeOutlet;
        bottomPresBCPS = &bottomAgeOutlet;
    }

	sim.rhoPSMean = rhoPSMean;

    slog(1, SLOG_INFO, "-> Setting up grid...\n");

	inletNodes = arraylist_new(L+1);
	inletCoords = arraylist_new(L+1);

	woundCoords = arraylist_new(50);	// sensible guess for wound size

	for(i=0; i<geometry->x; i++)
		for(j=0; j<geometry->y; j++) {
			int r,g,b;
			int idx = j*geometry->x+i;
			r = geometry->data[idx].r;
			g = geometry->data[idx].g;
			b = geometry->data[idx].b;
			
			Dynamics *dyn   = (Dynamics*)calloc(1, sizeof(Dynamics));
			Dynamics *dynPS = (Dynamics*)calloc(1, sizeof(Dynamics));

			// FLUID 
			if(r==255 && g==255 && b==255) {
				// Copy omega values for parallel execution
				fluidData *fd = (fluidData *)calloc(1, sizeof(fluidData));

				fd->omega = sim.omega;

				// if we would like non-Newtonian dynamics
				if(useNonNewt) {
					fd->useNonNewt = useNonNewt;
					fd->np = &sim.nonNewt;
				}

				fd->computeStressTensor=computeStressTensor;
                fd->useSmagorinsky=useSmagorinsky;
                fd->smagorinskyC=smagorinskyC;

				// dynamics of fluid
				dyn->dynamicsFun = dynamicsFun;

				dyn->selfData = (void*)fd;

                if(doPlatelets) {
                    plateletData *pp = (plateletData *) calloc(1, sizeof(plateletData));
                    pp->omega = sim.coupledOmega;
                    pp->margForceRatio = margForceRatio;
                    dynPS->selfData = (void *) pp;
                }
                else{
                    ageData *ad = (ageData*)calloc(1, sizeof(ageData));
                    ad->dt = sim.dt;
                    ad->fluidPop = sim.latticeA[i+1][j+1].fPop;
                    dynPS->selfData = (void *) ad;
                }

                // dynamics of coupled platelet field
				dynPS->dynamicsFun = dynamicsFunPS;

				
				// set geometry type descriptor
				setGeometry(&sim, i+1,j+1, FLUID);

			} // PFLUID
			else if(r==0 && g==255 && b==255) {
				// Copy omega values for parallel execution
				fluidData *fd = (fluidData *)calloc(1, sizeof(fluidData));

				fd->omega = sim.omega;

				// if we would like non-Newtonian dynamics
				if(useNonNewt) {
					fd->useNonNewt = useNonNewt;
					fd->np = &sim.nonNewt;
				}

				fd->computeStressTensor=computeStressTensor;
                fd->useSmagorinsky=useSmagorinsky;
                fd->smagorinskyC=smagorinskyC;

				// dynamics of fluid
				dyn->dynamicsFun = dynamicsFun;

				dyn->selfData = (void*)fd;

                if(doPlatelets) {
                    plateletData *pp = (plateletData *) calloc(1, sizeof(plateletData));
                    pp->omega = sim.coupledOmega;
                    // Protected fluid does not experience margination
                    pp->margForceRatio = 0;
                    dynPS->selfData = (void *) pp;
                }
                else{
                    ageData *ad = (ageData*)calloc(1, sizeof(ageData));
                    ad->dt = sim.dt;
                    ad->fluidPop = sim.latticeA[i+1][j+1].fPop;
                    dynPS->selfData = (void *) ad;
                }

                // dynamics of coupled platelet field at protected sites
				dynPS->dynamicsFun = dynamicsFunPSProt;


				// set geometry type descriptor
				setGeometry(&sim, i+1,j+1, PFLUID);
			}// NO-SLIP WALL
			else if(r==0 && g==0 && b==0) {

				dyn->dynamicsFun = &bounceBack;
				dynPS->dynamicsFun = &bounceBack;

				// set geometry type descriptor
				setGeometry(&sim, i+1,j+1, WALL);

			} // Velocity boundary
			else if(r==255 && g==0 && b==0) {

				VelocityBCData *velBC = (VelocityBCData *)calloc(1, sizeof(VelocityBCData));
				PressureBCData *presPSBC = (PressureBCData *)calloc(1, sizeof(PressureBCData));

				// Check side if this is the first velocity node
				//if(inletSide==UNSET) {	//commented out to prepare to have more inlets
					if(i == geometry->x-1) {
						inletSide = RIGHT;
						dyn->dynamicsFun = rightVelBC;
                        dynPS->dynamicsFun = rightVelBCPS;
					}
					else if (j == 0) {
						inletSide = TOP;
						dyn->dynamicsFun = topVelBC;
                        dynPS->dynamicsFun = topVelBCPS;
					}
					else if (i == 0) {
						inletSide = LEFT;
						dyn->dynamicsFun = leftVelBC;
                        dynPS->dynamicsFun = leftVelBCPS;
					}
					else if (j == geometry->y-1) {
						inletSide = BOTTOM;
						dyn->dynamicsFun = bottomVelBC;
                        dynPS->dynamicsFun = bottomVelBCPS;
					}
					else {
						slog(0, SLOG_WARN, "Wrong velocity boundary position at: (%d,%d)!!!\n", i,j);
						freeSimulation();
						exit(-1);
						// TODO: Should die gracefully
					}
				//}

				// Set up bulk dynamics for collision after the filling of unknown distribution elements
				fluidData *fd = (fluidData *)calloc(1, sizeof(fluidData));

				fd->omega = sim.omega;

				// if we would like non-Newtonian dynamics
				if(useNonNewt) {
					fd->useNonNewt = useNonNewt;
					fd->np = &sim.nonNewt;
				}

				fd->computeStressTensor=computeStressTensor;
                fd->useSmagorinsky=useSmagorinsky;
                fd->smagorinskyC=smagorinskyC;

				Dynamics *bulk = (Dynamics*)calloc(1, sizeof(Dynamics));

				bulk->dynamicsFun = dynamicsFun;

				bulk->selfData = (void*)fd;

				Dynamics *bulkPS = (Dynamics*)calloc(1, sizeof(Dynamics));

                    // Set BC properties
				velBC->bulkDynamics = bulk;
				velBC->ux = 0; velBC->uy = 0;	// Initialize these, just in case...

				// Set boundary condition info
				dyn->selfData = (void*)velBC;

                if(doPlatelets) {
                    plateletData *pp = (plateletData *) calloc(1, sizeof(plateletData));
                    pp->omega = sim.coupledOmega;
                    pp->margForceRatio = margForceRatio;

                    bulkPS->dynamicsFun = dynamicsFunPS;
                    bulkPS->selfData = (void *) pp;

                    presPSBC->bulkDynamics = bulkPS;
                    presPSBC->uPar = 0;
                    presPSBC->rho = rhoPSMean;

                    dynPS->selfData = (void *) presPSBC;
                }
                else{
                    ageData *ad = (ageData*)calloc(1, sizeof(ageData));
                    ad->dt = sim.dt;
                    ad->fluidPop = sim.latticeA[i+1][j+1].fPop;
                    dynPS->selfData = (void *) ad;
                }

				// store inlet node dynamics for time dependent settings
				arraylist_append(inletNodes, (void*)dyn);
				coord *c = (coord*)calloc(1, sizeof(coord));
				c->x=i+1; c->y=j+1;
				arraylist_append(inletCoords, (void*)c);

				// set geometry type descriptor
				setGeometry(&sim, i+1,j+1, BDFLUID);

			} // Pressure boundary
			else if(r==0 && g==255 && b==0) {
				
				PressureBCData *presBC = (PressureBCData *)calloc(1, sizeof(PressureBCData));
				PressureBCData *presPSBC = (PressureBCData *)calloc(1, sizeof(PressureBCData));

				if(i == geometry->x-1) {
					dyn->dynamicsFun = rightPresBC;
					dynPS->dynamicsFun = rightPresBCPS;
				}
				else if (j == 0) {
					dyn->dynamicsFun = topPresBC;
					dynPS->dynamicsFun = topPresBCPS;
				}
				else if (i == 0) {
					dyn->dynamicsFun = leftPresBC;
					dynPS->dynamicsFun = leftPresBCPS;
				}
				else if (j == geometry->y-1) {
					dyn->dynamicsFun = bottomPresBC;
					dynPS->dynamicsFun = bottomPresBCPS;
				}
				else {
					slog(0, SLOG_WARN, "Wrong pressure boundary position at: (%d,%d)!!!\n", i,j);
					freeSimulation();
					exit(-1);
					// TODO: Should die more gracefully
				}

				//printf("Pressure boundary position at: (%d,%d)!!!\n", i,j);

				// Set up bulk dynamics for collision after the filling of unknown distribution elements
				fluidData *fd = (fluidData *)calloc(1, sizeof(fluidData));

				fd->omega = sim.omega;

				// if we would like non-Newtonian dynamics
				if(useNonNewt) {
					fd->useNonNewt = useNonNewt;
					fd->np = &sim.nonNewt;
				}

				fd->computeStressTensor=computeStressTensor;
                fd->useSmagorinsky=useSmagorinsky;
                fd->smagorinskyC=smagorinskyC;

				Dynamics *bulk = (Dynamics*)calloc(1, sizeof(Dynamics));
				bulk->dynamicsFun = dynamicsFun;
				bulk->selfData = (void*)fd;

				// Set BC properties
				presBC->bulkDynamics = bulk;
				presBC->rho = defaultDensity*outletDefaultDensity; presBC->uPar = 0;	// Initialize these, just in case...

				// Set boundary condition info
				dyn->selfData = (void*)presBC;

                if(doPlatelets) {
                    plateletData *pp = (plateletData *) calloc(1, sizeof(plateletData));
                    pp->omega = sim.coupledOmega;
                    pp->margForceRatio = margForceRatio;

                    Dynamics *bulkPS = (Dynamics *) calloc(1, sizeof(Dynamics));
                    bulkPS->dynamicsFun = dynamicsFunPS;
                    bulkPS->selfData = (void *) pp;

                    presPSBC->bulkDynamics = bulkPS;
                    presPSBC->uPar = 0;
                    presPSBC->rho = rhoPSMean;

                    dynPS->selfData = (void *) presPSBC;
                }
                else {
                    ageData *ad = (ageData*)calloc(1, sizeof(ageData));
                    ad->dt = sim.dt;
                    ad->fluidPop = sim.latticeA[i+1][j+1].fPop;
                    dynPS->selfData = (void *) ad;
                }


                // set geometry type descriptor
				setGeometry(&sim, i+1,j+1, BDFLUID);

			} // Porous material
			else  if(r==0 && g==0 && b==255) {
				// Copy omega and porosity values for parallel execution
				porousData *pd = (porousData*)calloc(1, sizeof(porousData));

				pd->fp.omega = sim.omega;
				pd->c0 = c0;
				pd->c1 = c1;

				if(useNonNewt) {
					pd->fp.useNonNewt = useNonNewt;
					pd->fp.np = &sim.nonNewt;
				}

				pd->fp.computeStressTensor=computeStressTensor;
                pd->fp.useSmagorinsky=useSmagorinsky;
                pd->fp.smagorinskyC=smagorinskyC;

				// dynamics of fluid
                if(useMRT)
				    dyn->dynamicsFun = &mrtForced;
                else
                    dyn->dynamicsFun = &bgkForced;

				dyn->selfData = (void*)pd;

                if(doPlatelets){
                    plateletData *pp = (plateletData *)calloc(1,sizeof(plateletData));
                    pp->omega = sim.coupledOmega;
                    pp->margForceRatio = margForceRatio;

                    dynPS->selfData = (void*)pp;
                }
                else {
                    ageData *ad = (ageData*)calloc(1, sizeof(ageData));
                    ad->dt = sim.dt;
                    ad->fluidPop = sim.latticeA[i+1][j+1].fPop;
                    dynPS->selfData = (void *) ad;
                }

				// dynamics of coupled platelet field
				dynPS->dynamicsFun = dynamicsFunPS;

				// set geometry type descriptor
				setGeometry(&sim, i+1,j+1, FLUID);

			} // Wound location
			else if(r==255 && g==255 && b==0) {

				coord *c = (coord*)calloc(1, sizeof(coord));
				c->x=i+1; c->y=j+1;
				arraylist_append(woundCoords, (void*)c);

				// Otherwise treat it as a fluid cell

				// Copy omega values for parallel execution
				fluidData *fd = (fluidData *)calloc(1, sizeof(fluidData));

				fd->omega = sim.omega;

				// if we would like non-Newtonian dynamics
				if(useNonNewt) {
					fd->useNonNewt = useNonNewt;
					fd->np = &sim.nonNewt;
				}

				fd->computeStressTensor=computeStressTensor;
                fd->useSmagorinsky=useSmagorinsky;
                fd->smagorinskyC=smagorinskyC;

				// dynamics of fluid
				dyn->dynamicsFun = dynamicsFun;

				dyn->selfData = (void*)fd;

                if(doPlatelets){
                    plateletData *pp = (plateletData *)calloc(1,sizeof(plateletData));
                    pp->omega = sim.coupledOmega;
                    pp->margForceRatio = margForceRatio;
                    dynPS->selfData = (void*)pp;
                }
                else{
                    ageData *ad = (ageData*)calloc(1, sizeof(ageData));
                    ad->dt = sim.dt;
                    ad->fluidPop = sim.latticeA[i+1][j+1].fPop;
                    dynPS->selfData = (void *) ad;
                }

                // dynamics of the coupled field
                dynPS->dynamicsFun = dynamicsFunPS;

                // set geometry type descriptor
				setGeometry(&sim, i+1,j+1, FLUID);
			} // Error handling
			else {
				slog(0, SLOG_WARN, "Unsupported color caught in geometry at: (%d,%d)!!!!\n", i, j);
				slog(0, SLOG_WARN, "Undefined color: r:%d g:%d b:%d\n", r, g, b);
				slog(1, SLOG_INFO, "Exiting...\n");
				freeSimulation();
				exit(-1);
			}
			
			// Setting the dynamics
			setDynamics(&sim, i+1, j+1, dyn);
			setPSDynamics(&sim, i+1, j+1, dynPS);

			// Configuring the distribution functions to stationary initial state
			iniEquilibrium(&sim, i+1, j+1, defaultDensity, 0., 0.);
			iniPSEquilibrium(&sim, i+1, j+1, rhoPSMean, 0., 0.);

		}
	slog(1, SLOG_INFO, "-> Unused cells post processing...\n");

	int numUnused=0;
	for(i=1; i<=geometry->x; i++)
		for(j=1; j<=geometry->y; j++) {
			if(getGeometry(&sim, i,j) == WALL)
			{
				int k, l, neigh;
				neigh = 0;
				for(k=-1; k<=1;k++)
					for(l=-1;l<=1; l++){
						int cx, cy;
						cx = i+k;
						cy = j+l;
						if(cx > 0 && cx <= geometry->x && cy > 0 && cy <= geometry->y)
						{
							if(getGeometry(&sim, cx, cy)!= WALL && getGeometry(&sim, cx, cy)!= UNUSED)
								neigh++;
						}

				}
				if(neigh==0)
				{
					getDynamics(&sim.latticeA[i][j])->dynamicsFun=&nodynamics;
					getDynamics(&sim.psLatticeA[i][j])->dynamicsFun=&nodynamics;
					if(sim.doubleBuffered){
						getDynamics(&sim.latticeB[i][j])->dynamicsFun=&nodynamics;
						getDynamics(&sim.psLatticeB[i][j])->dynamicsFun=&nodynamics;
					}
					setGeometry(&sim, i,j, UNUSED);
					numUnused++;
				}
			}
		}

	slog(1, SLOG_INFO, "-> Found %d unused cells!\n", numUnused);

	slog(1, SLOG_INFO, "-> Boundary condition post processing...\n");

	real scale = getScale(inletBC, 0);

	setVelocityBoundary(inletNodes, L, inletSide, lbm_uMax *scale);

	// Reinitialize velocity boundary equilibrium to given BC data
	int m;
	for(m = 0; m < L; m++){
		Dynamics *d = (Dynamics *)inletNodes->data[m];
		VelocityBCData *bc = (VelocityBCData *)d->selfData;
		coord *c = (coord *)inletCoords->data[m];

		iniEquilibrium(&sim, c->x, c->y, 1., bc->ux, bc->uy);
		//TODO: Ini PS boundary too to the same velocity?
	}

	// Update global phase
	phase = 2;

	slog(1, SLOG_INFO, "-> Construction done successfully!\n");
	
	return 0;
}

int freeSimulation()
{
	slog(1, SLOG_INFO, "-> Cleaning up...\n");

	// Check if we are through initialization.
	if(phase > 0) {
		if(isCoupled==1 && calcThrombus == 1)
			freeNearWallCells(nearWallCells);

		arraylist_free(inletNodes);
		destructSim(&sim);
	}

	//TODO: Loop through geometry to deallocate all the dynamics!

	return 0;
}

void printCurrentInfo(Simulation *sim, const char *file, int iT, int iTmax, real uMaxPhys, const char *note)
{
    char model[60];

    if(useMRT)
        sprintf(model, "MRT");
    else
        sprintf(model, "BGK");

    if(useSmagorinsky)
        sprintf(model+3, " + Smagorinsky(%f)", smagorinskyC);

	FILE *f;

	f=fopen(file, "w");

	fprintf(f, "\n\t|##>    medFlow2D    <##|\n\n");

	fprintf(f,"Version: %d.%d\n", VERSION_MAJOR, VERSION_MINOR);
	fprintf(f,"Dynamics: %s\n", model);
	fprintf(f,"Geometry extents: x: %d, y: %d\n", sim->lx, sim->ly);
	fprintf(f,"Lattice sizes:\n\tdx: %e m\n\tdt: %e s\n\tdm: %e kg/m^3\n", sim->dx, sim->dt, sim->dm);
	fprintf(f,"Parameters:\n\tRe: %0.2f\n\t Omega: %f (nu: %f)\n", sim->Re, sim->omega, (1./3.)*(1./sim->omega-0.5));
	fprintf(f,"Maximal expected velocity: %f m/s\n", uMaxPhys);

	if(c0 != 0 || c1 != 0)
		fprintf(f, "Porous properties: c0: %e, c1: %e\n", c0, c1);

	if(sim->coupled)
		fprintf(f,"Coupled simulation.\n");

    if(doPlatelets)
        fprintf(f,"Platelet model is active. CoupledOmega: %f\n", sim->coupledOmega);

	if(sim->thrombus)
		fprintf(f, "Thrombus formation model is active.\n");

    if(doAgeing)
        fprintf(f,"Fluid aging calculation is active.\n");

	if(sim->useNonNewtonian)
		fprintf(f, "Non-Newtonian dynamics is enabled!\n");

	fprintf(f,"\nCurrent iteration: %d / %d\n", iT, iTmax);
	fprintf(f,"Current time: %f s / %f s\n", iT*sim->dt, iTmax*sim->dt);

	fprintf(f,"State: %s\n", note);

	fclose(f);
}

void saveCurrent(int iT, int iTMax, const char* note)
{
	saveVelPNG(&sim, filePath("currentVel.png", dirMonitor), lbm_uMax);
	saveRhoPNG(&sim, sim.activeLattice, filePath("currentRho.png", dirMonitor), imgRhoScaleMin, imgRhoScaleMax);		//scale: 50 - just small deviations
	saveRhoPNG(&sim, sim.activePsLattice, filePath("currentCoupledRho.png", dirMonitor), imgPSScaleMin, imgPSScaleMax);	//scale: 5 - large deviations by default
	printCurrentInfo(&sim, filePath("currentInfo.txt", dirMonitor), iT, iTMax, phys_uMax, note);
	saveFullSparseZ(&sim, filePath("currentVel.txt", dirMonitor), saveSparseStep);
}

// Get the fluid up to speed
int warmUp()
{

	// Initial velocity
	real startVel = getScale(inletBC, 0.0)* lbm_uMax;
	real maxULBMChangePerIt = sim.dt*maxULBMchangePerSec;
	int warmUpItNum = round(startVel / maxULBMChangePerIt);
	real warmUpChange = startVel/warmUpItNum;

	slog(1, SLOG_INFO, "-> %d iterations suggested for warmup!\n", warmUpItNum);

	int totCells = geometry->x*geometry->y;

	printf("\n\n");

	int iT;
	for (iT=0; iT<warmUpItNum; iT++) {

	        setVelocityBoundary(inletNodes, L, inletSide, iT*warmUpChange);

	        // on the right boundary, outlet condition grad_x u = 0
	        //updateZeroGradientBoundary();

	        // Update current image
	        if(iT% monitorUpdateInterval ==0) {

                loadBar(iT, warmUpItNum, 25, 25, totCells);

	        	saveCurrent(iT, warmUpItNum, "Warm up phase.");

                //check for external command
                switch(checkCommand(&sim, fullCmdFile))
                {
                    case -1:
                        slog(0, SLOG_INFO, "Remote EXIT signal received!          \n");

                        saveCurrent(iT, warmUpItNum, "Terminated by user.");

                        freeSimulation();

						slog(0, SLOG_INFO, "Exiting...\n");

                        exit(0);
                        break;

                    default:
                        break;

                }

	        	if(checkNanSparse(&sim, 10) < 0) {
                    //printf("\n\n!!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	        		slog(0, SLOG_ERROR, "Numerical simulation exploded at iteration: %d!\nCheck the ini parameters!\nQuitting...\n", iT);
	        		exit(-1);
	        	}
	        }

	        // Main fluid has to collideAA first, to update macroscopic velocity (NodeStore)
		    collideAA(&sim, sim.latticeA, numWorkers);

	        if(isCoupled) {
		        collideAA(&sim, sim.psLatticeA, numWorkers);  // Coupled passive scalar lattice
	        	propagateAA(&sim, sim.psLatticeA);
	        }

	        propagateAA(&sim, sim.latticeA);

	}

	printf("-> Warm up done.                                                 \n");
	return 0;
}

int simulate()
{
	int runIt = round(runTime/sim.dt);
	slog(1, SLOG_INFO, "-> Now running simulation for %f seconds (%d iterations)...\n", runTime, runIt);

	int saveIt = round(saveInterval /sim.dt);

	int totCells = geometry->x*geometry->y;

	printf("\n\n");

	int iT;
	for (iT=0; iT<runIt; iT++) {


		        real currentTime = iT*sim.dt;	// current physical time [s]
		        real currentUMax = getScale(inletBC, currentTime)* lbm_uMax;

		        setVelocityBoundary(inletNodes, L, inletSide, currentUMax);

		        // on the right boundary, outlet condition grad_x u = 0
		        //updateZeroGradientBoundary();

		        // Update current image
		        if(iT% monitorUpdateInterval ==0) {

                    loadBar(iT, runIt, 25, 25, totCells);

		        	saveCurrent(iT, runIt, "Simulation.");

		        	//check for external command
		        	switch(checkCommand(&sim, fullCmdFile))
                    {
                        case -1:
                            slog(0, SLOG_INFO, "Remote EXIT signal received!          \n");

                            saveCurrent(iT, runIt, "Terminated by user.");

                            freeSimulation();

                            slog(0, SLOG_INFO, "Exiting...\n");

                            exit(0);
                            break;

                        default:
                            break;

                    }

		        	if(checkNanSparse(&sim, 10) < 0) {
                        //printf("\n\n!!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                        slog(0, SLOG_ERROR, "Numerical simulation exploded at iteration: %d!\nCheck the ini parameters!\nQuitting...\n", iT);
		        		exit(-1);
		        	}
		        }

		        // Save current result
		        if ((iT-1)% saveIt ==0) {
                    char fName[80];

					if(doSaveImg)
					{
                        sprintf(fName, "%s/vel_%07d.png", dirResImg, iT);
                        saveVelPNG(&sim, fName, lbm_uMax);
                        //sprintf(fName, "%s/rho_%07d.ppm", dirResImg, iT);
                        //saveRhoPPM(&sim, sim.activeLattice, fName, 1.0, imgRhoScale);
                        sprintf(fName, "%s/coupledrho_%07d.png", dirResImg, iT);
                        saveRhoPNG(&sim, sim.activePsLattice, fName, imgPSScaleMin, imgPSScaleMax);
					}

		            if(doSaveData){
		            	sprintf(fName, "%s/full_%07d.txt", dirResData, iT);
		            	if(useCompression)
				            saveFullZ(&sim, fName);
			            else
				            saveFull(&sim, fName);
		            }

                    if(doSaveCoag && calcThrombus) {
                        sprintf(fName, "%s/coag_%07d.txt", dirResData, iT);
                        saveNearWallData(fName, nearWallCells);
                    }
		        }

		        // ******** Collisions ********

		        // Main fluid has to collideAA first, to update macroscopic velocity (NodeStore)
		        collideAA(&sim, sim.latticeA, numWorkers);

		        if(isCoupled) {
			        collideAA(&sim, sim.psLatticeA, numWorkers);  // Coupled passive scalar lattice

		        	// develop geometry
		        	if(calcThrombus)
		        		if (iT%coagulationUpdateInterval==0) {
		        			updateNearWallListHistory(&sim, sim.latticeA, nearWallCells);
		        			updateNearWallCells(&sim, nearWallCells);
		        		}
		        }

		        // ******** Propagations ********

		        propagateAA(&sim, sim.latticeA);

		        if(isCoupled)
		        	propagateAA(&sim, sim.psLatticeA);

		}

	saveCurrent(iT, runIt, "Simulation finished.");
	return 0;
}

int initThrombus()
{
	// Initialization of thrombus model

	slog(1, SLOG_INFO, "-> Initializing thrombus model...\n");

	nearWallCells = buildNearWallCells(&sim, sim.latticeA);    // reasonable bet on length

	coagulationUpdateInterval = round(activationTime/sim.dt/(real)getMaxHistory());

	// Modify wound cells to have higher coag. prob.
	int m;
	for(m=0; m < woundCoords->length; m++) {
		coord *c = (coord *)(woundCoords->data[m]);
		setRhoADP(nearWallCells, c->x, c->y, getDefaultRhoADP() * woundMultiplier);
	}

	return 0;
}

/*
 * Main entry point.
 */
int main(int argc, char** argv) {

	if (argc < 2) {
		printf("Usage: %s setup.ini\n", argv[0]);
		exit(-1);
	}

	fileInputINI = basename(strdup(argv[1]));
	dirInput = dirname(strdup(argv[1]));

	init_slog(NULL, 0, LOG_LEVEL);

	slog(1, SLOG_INFO, "**********************************\n");
	slog(1, SLOG_INFO, "*           medFlow2D            *\n");
	slog(1, SLOG_INFO, "**********************************\n\n");

	slog(1, SLOG_INFO, "-> Setting input directory to: %s\n", dirInput);
	slog(1, SLOG_INFO, "-> Setting input file to: %s\n", fileInputINI);

	fullCmdFile = filePath(cmdFile, dirInput);

	slog(1, SLOG_INFO, "-> Command file is: %s\n", fullCmdFile);

	// Clean the cmd file
	fclose(fopen(fullCmdFile, "w"));

	//Loading data from configuration file
	if (loadSetupIniFile(argv[1]) != 0) {
		slog(0, SLOG_ERROR, "Error during data loading!\nExiting gracefully...\n");
		exit(-1);
	}

	slog(1, SLOG_INFO, "-> Start logging to file (%s)...\n", filePath(fileLog, dirBase));

	init_slog(filePath(fileLog, dirBase), 1, LOG_LEVEL);

	slog(0, SLOG_INFO, "***** medFlow2D version %d.%d *****\n", VERSION_MAJOR, VERSION_MINOR);

	slog(1, SLOG_INFO, "-> Archiving setup.ini..\n");

	if(copyFile(filePath(fileInputINI, dirInput), filePath("setup.ini", dirBase)) == -1)
	{
		slog(0, SLOG_ERROR, "Cannot archive setup.ini!\nExiting gracefully...\n");
		exit(-1);
	}

	//Load transient boundary data
	slog(1, SLOG_INFO, "-> Loading time dependent boundary values...\n");

	inletBC = loadTransientBC(filePath(inletScaleFileName, dirInput));

	//Constructing the simulation
	if (constructSimulation() != 0) {
		slog(0, SLOG_ERROR, "Error during constructing the simulation!\nExiting gracefully...\n");
		exit(-1);
	}

	printSimDetails(&sim);

	// Check and create directory structure

	if (dirExists(dirBase) != 1) {
		slog(0, SLOG_ERROR, "Valid directory required for output saving!\n (Wrong working folder: %s)\nExiting..", dirBase);
		exit(-1);
	}

	// Directories
	sprintf(dirMonitor, "%s/monitor", dirBase);
	sprintf(dirResData, "%s/results_data", dirBase);
	sprintf(dirResImg, "%s/results_img", dirBase);

	if (dirExists(dirMonitor) != 1) {
		if (createDir(dirMonitor) != 0) {
			slog(0, SLOG_ERROR, "Cannot create folder for monitoring info!\n (Wrong folder: %s)\nExiting..", dirMonitor);
			exit(-1);
		}
	}

	if (doSaveData > 0) {
		if (dirExists(dirResData) != 1) {
			if (createDir(dirResData) != 0) {
				slog(0, SLOG_ERROR, "Cannot create folder for transient result files!\n (Wrong folder: %s)\nExiting..", dirResData);
				exit(-1);
			}
		}
	}

	if (doSaveImg > 0) {
		if (dirExists(dirResImg) != 1) {
			if (createDir(dirResImg) != 0) {
				slog(0, SLOG_ERROR, "Cannot create folder for image snapshots!\n (Wrong folder: %s)\nExiting..", dirResImg);
				exit(-1);
			}
		}
	}


	printBC(inletBC);


	if(calcPorousConsts){
		slog(1, SLOG_INFO, "-> Calculating porous model parameters..\n");
		calcPorousConstants(&sim, phys_nu, 0.8, 0.001, &c0, &c1);
	}

	slog(1, SLOG_INFO, "Porous material details:\t c0: %e\t c1: %e\n", c0, c1);

	// Have a different particle dynamics next to walls (no forcing term).
	// Also, it can make the simulation more stable
	modifyNearWallCellsPS(&sim, glycocalyxWidth, &bgkPS);

	// Save initial state
	slog(0, SLOG_INFO, "*** Simulation setup is ready. Saving inital state...\n");
	saveFullZ(&sim, filePath("initial.txt", dirBase));
	saveVelPNG(&sim, filePath("initialVel.png", dirBase), lbm_uMax);
	saveRhoPNG(&sim, sim.activeLattice, filePath("initialRho.png", dirBase), imgRhoScaleMin, imgRhoScaleMax);
	saveRhoPNG(&sim, sim.activePsLattice, filePath("initialCoupledRho.png", dirBase), imgPSScaleMin, imgPSScaleMax);


	if(useCheckPoint) {
		slog(1, SLOG_INFO, "-> Checkpointing enabled.\n");

		if(fileExists(filePath(fileCheckpoint, dirBase)))
		{	// Load checkpoint
			slog(1, SLOG_INFO, "-> Loading existing checkpoint...\n");
			loadCheckPoint(&sim, filePath(fileCheckpoint, dirBase));
		}
		else
		{	// Save checkpoint
			slog(1, SLOG_INFO, "-> No existing checkpoint data found!\n");

			warmUp();

			slog(1, SLOG_INFO, "-> Saving checkpoint...\n");
			saveCheckPoint(&sim, filePath(fileCheckpoint, dirBase));
		}
	}
	else
		warmUp();

	// TODO: Probe resulting platelet density somewhere and set is as boundary condition?

	if(calcThrombus) {
		initThrombus();
	}

	simulate();

	freeSimulation();

	slog(1, SLOG_INFO, "-> Simulation done. Exiting.\n");

    return (0);
}


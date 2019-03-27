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

#ifndef LB_H
#define LB_H

//#include "stddef.h"
#include "../common.h"

// Should use single character numbers :)
#define VERSION_MAJOR 1
#define VERSION_MINOR 1    // or not...

/* D2Q9 lattice constants                                        */
/*****************************************************************/

  // lattice weights
static const real t[9] = { 4./9., 1./9., 1./36., 1./9., 1./36.,
                             1./9., 1./36., 1./9., 1./36. };
  // lattice velocities
static const int c[9][2] = {
    { 0, 0},
    { 1, 0}, { 1, -1}, {0, -1}, {-1, -1},
    {-1, 0}, {-1,  1}, {0,  1}, { 1,  1}
};

  // opposite directions, for bounce back implementation
static const int oppositeOf[9] = { 0, 5, 6, 7, 8, 1, 2, 3, 4 };

typedef enum {
	WALL,			// regular bounceback wall
	FLUID,			// regular fluid 
	BDFLUID,		// boundary fluid cells (ex.: ZouHe)
	NWFLUID,		// NearWall Fluid - can turn into wall
    COAGWALL,       // Wall created from coagulated fluid
	PFLUID,			// Protected fluid cell (no special process, like coagulation can happen here)
	UNUSED			// Unused wall cell with only wall neighbours
} NodeType;


/* global state modifiers                                       */
/****************************************************************/
//void setComputeStressTensor(bool value);
//inline bool getComputeStressTensor();

/* struct NodeStore												*/
/****************************************************************/
/* Easily extendable store structure to help reuse data. 
 * Detached to conserve memory. 
 */

typedef struct {
	real ux, uy;			    // Store macroscopic velocity
	real fx,fy;					// Store macroscopic force
    real Sij[3];              	// rate of strain tensor
    real mu;					// store dynamic viscosity
} NodeStore;


/* struct Dynamics                                               */
/*****************************************************************/
/*   emulation of a class that defines the dynamics of a lattice
 *   site. The function dynamicsFun contains the algorithm of the
 *   collision, and selfData points to the function arguments
 *   (for example the local viscosity)
 */
typedef struct {
    void   (*dynamicsFun) (real* fPop, void* selfData, NodeStore* store);
    void    * selfData;
} Dynamics;


/* struct Node                                                   */
/*****************************************************************/
/*  a lattice node, containing the data (distribution functions)
 *  for a D2Q9 lattice, and a pointer to the local dynamics, i.e.
 *  the collision term. Two "methods" are added to construct and
 *  initialize a node.
 */
typedef struct {
    real    fPop[9];
    Dynamics* dynamics;
} Node;

// Structures extending normal fluid properties
typedef struct {
	real nu_0;
	real nu_inf;
	real lambda;
	real a;
	real n;
} NonNewtParam;

typedef struct {
	real omega;
	int useNonNewt;
	bool computeStressTensor;
	bool useSmagorinsky;
	NonNewtParam *np;
    real smagorinskyC;
} fluidData;

typedef struct {
	real omega;
	real margForceRatio;
} plateletData;

typedef struct {
	real c0;
	real c1;
	fluidData fp;
} porousData;

typedef struct {
    real dt;
    real *fluidPop;
} ageData;

/* struct Simulation                                             */
/*****************************************************************/
/*  a full D2Q9 LB simulation, containing the allocated memory
 *  for the lattice. Some "methods" are added to initiate the
 *  dynamics.
 */
typedef struct {
	int 		 lx, ly;				// Simulation grid resolution
	int			 coupled;				// Use tmpLattice(0/1)?
	int 		 thrombus;				// Use thrombus model (0/1)?
	int 		 useNonNewtonian; 		// Use non-Newtonian dynamics
    int          useSmagorinsky;        // Use Smagorinsky turbulence modelling
    int          doubleBuffered;        // (0/1) for AA or AB type addressing
	real 		 dx, dt, dm;			// Scale: LBM values in SI metric
	real 		 Re;					// Reynolds
	real 	 	 omega;					// Relaxation frequency
	real	 	 coupledOmega;			// Coupled field relaxation frequency
	real		 rhoPSMean;				// mean concentration of passive scalar field
	NodeStore**  store;					// Macroscopic velocity
	NodeStore*	 storeMemoryChunk;		// contiguous raw memory for UV (probably rho later)
    Node*  		 memoryChunkA;       	// contiguous raw memory
    Node*        psMemoryChunkA;  		// contiguous raw tmp memory
	Node** 		 latticeA;           	// lattice, points to raw memory
	Node**       activeLattice;         // Which lattice to read from in case of double buffering
    Node**       psLatticeA;        	// tmp lasttice, points to raw memory (passiveScalar)
	Node**       activePsLattice;       // Which tmpLattice to read from in case of double buffering
    Node*  		 memoryChunkB;       	// contiguous raw memory
    Node*        psMemoryChunkB;  		// contiguous raw tmp memory
    Node** 		 latticeB;           	// lattice, points to raw memory
    Node**       psLatticeB;        	// tmp lasttice, points to raw memory (passiveScalar)
    NodeType*    geometryMemoryChunk;	// Geometry flag
    NodeType**	 geometry;				// Geometry flag
    NonNewtParam nonNewt;				// Non-Newtonian dynamics parameters
} Simulation;

void constructSim(Simulation* sim, int lx, int ly, int doDoubleBuffer, int coupled, int thrombus);  // allocate memory for a full simulation ("constructor")
void destructSim(Simulation* sim); // free the memory for the simulation ("destructor")

void constructNode(Node* node);
void iniEquilibrium(Simulation *sim, int x, int y, real rho, real ux, real uy);
void iniPSEquilibrium(Simulation *sim, int x, int y, real rho, real ux, real uy);

void flipActiveLattice(Simulation *sim);

/* struct for parallel execution                            */
/************************************************************
 * Every worker receives this data structure.
 */

typedef struct {
    Simulation* sim;        // simulation pointer
    Node**      nodes;      // Nodes to compute in simulation
    int         id;         // Worker ID
	int 		fromIdx;
	int 		toIdx;
} ParallelData;

void setMargForceRatio(real mFR);
real getMargForceRatio ();
void cpyDistribution(Simulation *sim, Node **latticeA, Node **latticeB);
void setDynamics(Simulation *sim, int x, int y, Dynamics *dyn);
void setPSDynamics(Simulation *sim, int x, int y, Dynamics *dyn);
Dynamics* getDynamics(Node *cell);
void setGeometry(Simulation *sim, int iX, int iY, NodeType type);
NodeType getGeometry(Simulation *sim, int iX, int iY);
void collideAA(Simulation *sim, Node **lattice, int numWorkers);
void propagateAA(Simulation *sim, Node **lattice);
void makePeriodic(Simulation* sim, Node** lattice);

void calculateShearStressMaxDomain(Simulation* sim, real** shearMax, real** shearAngle);
inline void calculateShearStressMax(NodeStore* store, real* shearMax, real* shearAngle);

/* some free helper functions                                    */
/*****************************************************************/

//Compute macroscopic density from population
void inline computeRho(real* f, real* rho);

inline real computePSEquilibrium(int iPop, real rho, real ux, real uy);

  // compute density and velocity from the f's 
inline void computeMacros(real* f, real* rho, real* ux, real* uy);

  // compute local equilibrium from rho and u
inline real computeEquilibrium(int iPop, real rho, real ux, real uy, real uSqr);

  // Passive scalar with simplified BGK collision
void bgkPS(real *fPop, void *selfData, NodeStore* store);

  // Passive scalar with simplified BGK collision plus a force term
void bgkForcedPS(real *fPop, void *selfData, NodeStore* store);

  // bgk collision term
void bgk(real* fPop, void* selfData, NodeStore* store);

  // mrt collision term
void mrt(real* fPop, void* selfData, NodeStore* store);

  // forced mrt modell for porous media
void  mrtForced(real *fPop, void *selfData, NodeStore *store);

  // equilibrium boundary condition dynamics
void equBC(real * fPop, void* selfData, NodeStore* store);

  // porous media collision
void bgkForced(real *fPop, void *selfData, NodeStore *store);

  // ageing dynamics
void ageingDyn(real *fPop, void *selfData, NodeStore *store);

#endif

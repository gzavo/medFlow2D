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

//#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <libgen.h>
#include "slog.h"
#include "utils.h"

// Include compression library
#include "miniz.h"
#include "../lbm/lb.h"

extern real inline getScale(transientBC *bc, real phys_time);

static inline real getJetR(real gray);
static inline real getJetG(real gray);
static inline real getJetB(real gray);

inline real interpolate1(real x1, real x2, real xi, real y1, real y2)
{
    real r = (xi-x1)/(x2-x1);
    return y1+r*(y2-y1);
}

void fprintSimDetails(Simulation *sim, char *fileName) {
    FILE *f;
    f = fopen(fileName, "w");

    fprintf(f, "\n************************************\n");
    fprintf(f, "*           medFlow2D %d.%d         *\n", VERSION_MAJOR, VERSION_MINOR);
    fprintf(f, "************************************\n\n");

    fprintf(f, "Domain -> x: %d\t y:%d\n", sim->lx, sim->ly);
    printf("Lattice units -> dx: %e\t dt: %e\t dm: %e\n", sim->dx, sim->dt, sim->dm);
    fprintf(f, "Parameters -> Re: %0.2f\t Omega: %f (nu: %f)\n", sim->Re, sim->omega, (1. / 3.) * (1. / sim->omega - 0.5));

    if(sim->useSmagorinsky)
        fprintf(f, "Smagorinsky turbulence model is enabled!\n");

    if (sim->coupled)
        fprintf(f, "Coupled simulation. CoupledOmega: %f\n\n", sim->coupledOmega);

    if (sim->useNonNewtonian)
        fprintf(f, "Non-Newtonian fluid material is enabled!\n");

    fclose(f);
}

void printSimDetails(Simulation *sim) {

    slog(1, SLOG_INFO, "*********  Parameters  ************\n");
    slog(1, SLOG_INFO, "Domain -> x: %d\t y:%d\n", sim->lx, sim->ly);
    slog(1, SLOG_INFO, "Lattice units -> dx: %e\t dt: %e\t dm: %e\n", sim->dx, sim->dt, sim->dm);
    slog(1, SLOG_INFO, "Parameters -> Re: %0.2f\t Omega: %f (nu: %f)\n", sim->Re, sim->omega, (1. / 3.) * (1. / sim->omega - 0.5));

    if(sim->useSmagorinsky)
        slog(1, SLOG_INFO, "Smagorinsky turbulence model is enabled!\n");

    if (sim->coupled)
        slog(1, SLOG_INFO, "Coupled simulation. CoupledOmega: %f\n\n", sim->coupledOmega);

    if (sim->useNonNewtonian)
        slog(1, SLOG_INFO, "Non-Newtonian fluid material is enabled!\n");

    fflush(stdout);
}


int createDir(const char *path)
{

#ifdef _WIN32 // note the underscore: without it, it's not msdn official!
    // Windows (x64 and x86)
    if(mkdir(path) != 0)
        return 1;
#elif __unix__ // all unices, not all compilers
    // Unix
    if(mkdir(path, 0777) != 0)
        return 1;
#elif __linux__
    // linux
    if(mkdir(path, 0777) != 0)
        return 1;
#elif __APPLE__
    // Mac OS, not sure if this is covered by __posix__ and/or __unix__ though...
    if(mkdir(path, 0777) != 0)
        return 1;
#endif

   return 0;
}

int dirExists(const char *path)
{
    struct stat info;

    if(stat( path, &info ) != 0)
        return 0;
    else if(info.st_mode & S_IFDIR)
        return 1;
    else
        return 0;
}

// Check if file available for reading (might not for writing tough)
int fileExists(const char *fname)
{
    FILE *file= fopen(fname, "r");
    if (file)
    {
        fclose(file);
        return 1;
    }
    return 0;
}

// base dir extrection
char *baseName(char const *path)
{
    char *s = strrchr(path, '/');

    if (!s)
    	return strdup(path);
    else
        return strdup(s + 1);
}

int copyFile(char *old_filename, char  *new_filename)
{
    FILE  *ptr_old, *ptr_new;

    ptr_old = fopen(old_filename, "r");

    if(!ptr_old)
        return  -1;

    ptr_new = fopen(new_filename, "w");

    if(!ptr_new)
    {
        fclose(ptr_old);
        return  -1;
    }

    int ch;
    while( ( ch = fgetc(ptr_old) ) != EOF )
        fputc(ch, ptr_new);

    fclose(ptr_new);
    fclose(ptr_old);
    return  0;
}

// check command file
int checkCommand(Simulation *sim, const char * fileName)
{
	char signal[80];
	char cmdStr[256];

//	Command *cmd = (Command*)calloc(1, sizeof(Command));
//	cmd->commandStr = (char*)calloc(256, sizeof(char));

	if(!fileExists(fileName))
		return 0;

	struct stat st;
	stat(fileName, &st);
	int size = st.st_size;

	if(size > 0) {
		FILE *file;
		file = fopen(fileName, "r");
		fscanf(file, "%s %s", signal, cmdStr);
		fclose(file);

		//printf("Command found %s, %s\n", signal, cmdStr);

		if(strncmp("NOCMD", signal, 5)==0)
		{ /* Skip */ }
		else if (strncmp("SAVEF", signal, 5)==0)
    	{  /* TODO: write save F function */ }
		else if (strncmp("SAVEV", signal, 5)==0)
			saveVel(sim, cmdStr);
		else if (strncmp("SAVED", signal, 5)==0)
			saveRho(sim, sim->activeLattice, cmdStr);
		else if (strncmp("SAVEA", signal, 5)==0)
			saveFull(sim, cmdStr);
        else if (strncmp("EXIT", signal, 4)==0) {
            fclose(fopen(fileName, "w"));
            return -1;      // Exit command
        }

		// Command processed, erase file
		fclose(fopen(fileName, "w"));

		return 1;           // There was some command
	}

	return 0;               // There was no command
}

/*
 * Compute Poiseuille profile
 */
real poiseuilleProfile(int r, real uMax, int R)
{
	return uMax*(1.0f-( (real)(r*r)/ (real)(R*R)));
}

// Read in transient boundary condition file and fill structure
transientBC *loadTransientBC(const char *fileName) {

	transientBC *bc = (transientBC *)calloc(1, sizeof(transientBC));

    FILE *inletf;

    inletf = fopen(fileName, "r");
    // TODO: check if opened successfully

    if (inletf == NULL) {
        slog(0, SLOG_ERROR, "Error reading transient boundary data file!\nExiting...\n");
        exit(-1);
    }

    if(fscanf(inletf, "%d\n", &bc->numVal)==EOF) {
        slog(0, SLOG_ERROR, "Unrecognised transient boundary file format!\nExiting...\n");
        exit(-1);
    }

    bc->scales = (real*)calloc(bc->numVal, sizeof(real));
    bc->times = (real*)calloc(bc->numVal, sizeof(real));

    int i;
    for (i=0; i< bc->numVal; i++) {
        if(fscanf(inletf, "%lf %lf\n", &bc->times[i], &bc->scales[i])==EOF) {
        	slog(0, SLOG_ERROR, "Error: not enough data points in transient boundary file (%s)!", fileName);
        	exit(-1);
        }

    }

    fclose(inletf);

    bc->time_max = bc->times[bc->numVal-1];
    bc->time_current = 0;
    bc->currentCycle = 0;

    return bc;
}

void printBC(transientBC *bc)
{
	slog(1, SLOG_INFO, "Transient boundary condition info:\n");
	slog(1, SLOG_INFO, "  Number of values: %d\n", bc->numVal);
	slog(1, SLOG_INFO, "  Time span: %lf\n", bc->time_max);
	slog(1, SLOG_INFO, "  Currently active: %d\n", bc->time_current);
	slog(1, SLOG_INFO, "  Current cycle: %d\n", bc->currentCycle);
}

// Check for passive scalar field validity at chosen points
int checkPSSanity(Simulation *sim)
{
	// These two points seem OK
    int x1 = sim->lx/3;
    int x2 = 2*x1;
    int y = sim->ly/2;

    if(isnan(sim->activePsLattice[x1][y].fPop[0]) || isnan(sim->activePsLattice[x2][y].fPop[0]))
    {
        //Fck, at least exit gracefully.
        slog(0, SLOG_ERROR, "Passive scalar field diverged, should exit now...\n");
        return -1;
    }

    return 0;
}

// linear interpolation between distinct boundary values
real inline getScale(transientBC *bc, real phys_time)
{
	int quo = floor(phys_time/bc->time_max);
    real rem = fmod(phys_time, bc->time_max);

    if (quo > bc->currentCycle) {
    	bc->time_current = 0;
    	bc->currentCycle = quo;
    }
    else {
   		if (rem > bc->times[bc->time_current+1]) {
   			bc->time_current++;
   			if (bc->time_current >= (bc->numVal-1))
   				bc->time_current = 0;
    	}
    }

    real r = (rem-bc->times[bc->time_current])/(bc->times[bc->time_current+1]-bc->times[bc->time_current]);
    real intp = bc->scales[bc->time_current]+r*(bc->scales[bc->time_current+1]-bc->scales[bc->time_current]);
    //real intp = interpolate1(inlet_time[inlet_current], inlet_time[inlet_current+1], ctime, inlet_scale[inlet_current], inlet_scale[inlet_current+1]);
    //printf("\n %lf %lf %lf %lf", inlet_time[inlet_current], inlet_time[inlet_current+1], ctime, interpolate1(inlet_time[inlet_current], inlet_time[inlet_current+1], ctime, inlet_scale[inlet_current], inlet_scale[inlet_current+1]));
    return intp;
}


  // Progress bar code
void loadBar(int x, int n, int r, int w, int cells)
{
    // Calculate MLUPS and ETA
    static int lastIt = 0;
    static time_t lastTime = 0;

    time_t cTime = time(NULL);

    float dt = cTime - lastTime;

    if(dt <= 1.0 ) return;

    lastTime = cTime;
    // Only update r times.
    //if ( x % (n/r) != 0 ) return;

    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;

    float lups = ((x - lastIt)/dt)*cells;
    lastIt = x;
    long eta;
    if(lups == 0.0)
        eta = 0;
    else
        eta = ((n-x)*(cells/1000))/(lups/1000);

    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    int i;
    for (i=0; i<c; i++)
       printf("=");
 
    for (i=c; i<w; i++)
       printf(" ");
 
    printf("]");

    int rem1;
    int hours = eta / 3600;
    rem1 = eta % 3600;
    int min = rem1 / 60;
    int sec = rem1 % 60;

    printf(" MLUPS: %2.2f, ETA: %02d:%02d:%02d        ", lups/1e6, hours, min, sec);

    if(ratio <1.0)
        printf("\r");

    fflush(stdout);
}


// Various data saving routines

void saveCheckPoint(Simulation* sim, const char *fName) {
    FILE* oFile = fopen(fName, "wb");

    int iX, iY, iP;
    for (iY=0; iY<=sim->ly; iY++)
       for (iX=0; iX<=sim->lx; iX++)
    	   fwrite(sim->activeLattice[iX][iY].fPop, sizeof(real), 9, oFile);
//    	   for(iP=0; iP<9; iP++)
//    		   fprintf(oFile, "%f ", sim->lattice[iX][iY].fPop[iP]);

    for (iY=0; iY<=sim->ly; iY++)
       for (iX=0; iX<=sim->lx; iX++)
    	   fwrite(sim->activePsLattice[iX][iY].fPop, sizeof(real), 9, oFile);
//    	   for(iP=0; iP<9; iP++)
//    		   fprintf(oFile, "%f ", sim->tmpLattice[iX][iY].fPop[iP]);

    fclose(oFile);
}

void loadCheckPoint(Simulation* sim, const char *fName) {
    FILE* oFile = fopen(fName, "rb");

    int iX, iY, iP;
    for (iY=0; iY<=sim->ly; iY++)
       for (iX=0; iX<=sim->lx; iX++)
    	   fread(sim->activeLattice[iX][iY].fPop, sizeof(real), 9, oFile);
//    	   for(iP=0; iP<9; iP++)
//    		   fscanf(oFile, "%f ", &sim->lattice[iX][iY].fPop[iP]);

    for (iY=0; iY<=sim->ly; iY++)
       for (iX=0; iX<=sim->lx; iX++)
    	   fread(sim->activePsLattice[iX][iY].fPop, sizeof(real), 9, oFile);
//    	   for(iP=0; iP<9; iP++)
//    		   fscanf(oFile, "%f ", &sim->tmpLattice[iX][iY].fPop[iP]);

    fclose(oFile);
}


  // save the velocity field (norm) to disk
void saveVel(Simulation* sim, char fName[]) {
    FILE* oFile = fopen(fName, "w");
    int iX, iY;
    real ux, uy, uNorm, rho;
    for (iY=1; iY<=sim->ly; iY++) {
       for (iX=1; iX<=sim->lx; iX++) {
   		   computeMacros(sim->activeLattice[iX][iY].fPop, &rho, &ux, &uy);
   		   uNorm = sqrt(ux*ux+uy*uy);
   		   fprintf(oFile, "%f ", uNorm);
       }
       fprintf(oFile, "\n");
    }
    fclose(oFile);
}

// Helper function to avoid some chrash on error
int bytes_added( int result_of_sprintf )
{
    return (result_of_sprintf > 0) ? result_of_sprintf : 0;
}

// save the velocity field (sparsely to disk)
void saveFullSparse(Simulation* sim, char fName[], int step) {
  FILE* oFile = fopen(fName, "w");
  int iX, iY;
  int maxX = floor(sim->lx/step);
  int maxY = floor(sim->ly/step);

  fprintf(oFile, "%d\t%d\n", maxX, maxY);

  real ux, uy, fx, fy, rho, c, p;
  for (iY=0; iY<maxY; iY++) {
     for (iX=0; iX<maxX; iX++) {
    	 // Check if it is a wall node or not
    	 if(sim->geometry[iX*step+1][iY*step+1]!=WALL && sim->geometry[iX*step+1][iY*step+1] != COAGWALL) {
    		 //computeMacros(sim->lattice[iX*step+1][iY*step+1].fPop, &rho, &ux, &uy);
    		 computeRho(sim->activeLattice[iX*step+1][iY*step+1].fPop, &rho);
    		 computeRho(sim->activePsLattice[iX*step+1][iY*step+1].fPop, &c);
    		 //p = (rho-1.0)/3.0*sim->dm/sim->dx/sim->dt/sim->dt;
             //p = (rho-1.0)/3.0*sim->dx*sim->dx/sim->dt/sim->dt;
             p = (rho-1.0)/3.0*sim->dm/sim->dt/sim->dt;
    		 //rho = rho *sim->dm/sim->dx/sim->dx/sim->dx;
    		 ux = sim->store[iX*step+1][iY*step+1].ux *sim->dx/sim->dt;
    		 uy = sim->store[iX*step+1][iY*step+1].uy *sim->dx/sim->dt;
    		 fx = sim->store[iX*step+1][iY*step+1].fx *sim->dm*sim->dx/sim->dt/sim->dt;
    		 fy = sim->store[iX*step+1][iY*step+1].fy *sim->dm*sim->dx/sim->dt/sim->dt;

    		 //uNorm = sqrt(ux*ux+uy*uy);

    		 fprintf(oFile, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", ux, uy, rho, p, c-sim->rhoPSMean, fx, fy);
    	 }
    	 else	//If it is a wall cancel out all values
    		 fprintf(oFile, "%f\t%f\t%f\t%f\t%f\t%e\t%e\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
     }
     //fprintf(oFile, "\n");
  }
  fclose(oFile);
}

void saveFullSparseZ(Simulation* sim, char fName[], int step)
{

    // Variables
    char *strOut;
    int MAX_ELEMENT_LENGTH = 20;
    int MAX_ELEMENTS = 7;
    int length = 0;
    char *comment = "Output produced by medFlow2D";

    // Get file name parts
    char *baseFile = basename(strdup(fName));
    char *baseDir = dirname(strdup(fName));
    char archiveName[1024];
    char fullArchiveNAme[2048];

    // zip archive name
    sprintf(archiveName, "%s.zip", baseFile);
    sprintf(fullArchiveNAme, "%s/%s", baseDir, archiveName);

    // Sparse bounds
    int maxX = floor(sim->lx/step);
    int maxY = floor(sim->ly/step);

    //Guess size (over approximate)
    size_t length_approx = (size_t)(maxX*maxY*(MAX_ELEMENTS + MAX_ELEMENTS*MAX_ELEMENT_LENGTH));

    //slog(5, SLOG_DEBUG, "Allocation output buffer of length %d...\n", length_approx);

    //allocate buffer
    strOut = calloc(length_approx, sizeof(char));

    //slog(5, SLOG_DEBUG, "Filling up buffer (virtual %s)...\n", baseFile);

    // Write resolution
    length += bytes_added(sprintf(strOut+length, "%d\t%d\n", maxX+1, maxY+1));

    // Write output to buffer
    real ux, uy, fx, fy, rho, p, c;

    int iXF, iYF;
    for (iYF=0; iYF<=maxY; iYF++) {
        for (iXF=0; iXF<=maxX; iXF++) {
            int iX = iXF*step+1;
            int iY = iYF*step+1;

            if(sim->geometry[iX][iY]!=WALL && sim->geometry[iX][iY] != COAGWALL) {
                //computeMacros(sim->lattice[iX*step+1][iY*step+1].fPop, &rho, &ux, &uy);
                computeRho(sim->activeLattice[iX][iY].fPop, &rho);
                computeRho(sim->activePsLattice[iX][iY].fPop, &c);
                //p = ((rho-1.0)/3.0)*sim->dm/sim->dx/sim->dt/sim->dt;
                p = (rho - 1.0) / 3.0 * sim->dm / sim->dt / sim->dt;
                //p = (rho - 1.0) / 3.0 * sim->dx * sim->dx / sim->dt / sim->dt;
                //rho = rho *sim->dm/sim->dx/sim->dx/sim->dx;
                ux = sim->store[iX][iY].ux * sim->dx / sim->dt;
                uy = sim->store[iX][iY].uy * sim->dx / sim->dt;
                fx = sim->store[iX][iY].fx * sim->dm * sim->dx / sim->dt / sim->dt;
                fy = sim->store[iX][iY].fy * sim->dm * sim->dx / sim->dt / sim->dt;

                //uNorm = sqrt(ux*ux+uy*uy);
                length += bytes_added(
                        sprintf(strOut + length, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", ux, uy, rho, p, c - sim->rhoPSMean, fx,
                                fy));
                //fprintf(oFile, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", ux, uy, rho, sim->store[iX][iY].shear_max, sim->store[iX][iY].shear_ang, fx, fy);
            }
            else
                length += bytes_added(
                        sprintf(strOut + length, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        }
    }

    //slog(5, SLOG_DEBUG, "Beginning compression of output data (%d bytes) to (%s)\n", length, fullArchiveNAme);

    // If file exists, remove it
    if(fileExists(fullArchiveNAme))
        remove(fullArchiveNAme);

    // Compression
    mz_bool status;
    status = mz_zip_add_mem_to_archive_file_in_place(fullArchiveNAme, baseFile, strOut, length + 1, comment, (short)strlen(comment), MZ_BEST_COMPRESSION);
    if (!status)
    {
        slog(0, SLOG_ERROR, "Failed to compress output file!\n");
        //return EXIT_FAILURE;
    }

    //slog(5, SLOG_DEBUG, "Compression is done.\n");
    // Free up string buffer
    free(strOut);
}

void saveFull(Simulation* sim, char fName[]) {
  FILE* oFile = fopen(fName, "w");
  int iX, iY;

  fprintf(oFile, "%d\t%d\n", sim->lx, sim->ly);

  real ux, uy, fx, fy, rho, p, c;
  for (iY=1; iY<=sim->ly; iY++) {
     for (iX=1; iX<=sim->lx; iX++) {
   		 //computeMacros(sim->lattice[iX*step+1][iY*step+1].fPop, &rho, &ux, &uy);
   		 computeRho(sim->activeLattice[iX][iY].fPop, &rho);
   		 computeRho(sim->activePsLattice[iX][iY].fPop, &c);
   		 //p = ((rho-1.0)/3.0)*sim->dm/sim->dx/sim->dt/sim->dt;
         p = (rho-1.0)/3.0*sim->dm/sim->dt/sim->dt;
         //p = (rho - 1.0) / 3.0 * sim->dx * sim->dx / sim->dt / sim->dt;
   		 //rho = rho *sim->dm/sim->dx/sim->dx/sim->dx;
   		 ux = sim->store[iX][iY].ux *sim->dx/sim->dt;
   		 uy = sim->store[iX][iY].uy *sim->dx/sim->dt;
   		 fx = sim->store[iX][iY].fx *sim->dm*sim->dx/sim->dt/sim->dt;
   		 fy = sim->store[iX][iY].fy *sim->dm*sim->dx/sim->dt/sim->dt;

   		 //uNorm = sqrt(ux*ux+uy*uy);
   		 fprintf(oFile, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", ux, uy, rho, p, c-sim->rhoPSMean, fx, fy);
   		//fprintf(oFile, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", ux, uy, rho, sim->store[iX][iY].shear_max, sim->store[iX][iY].shear_ang, fx, fy);
     }
     //fprintf(oFile, "\n");
  }
  fclose(oFile);
}



// Compressed output
void saveFullZ(Simulation* sim, char fName[]) {

    // Variables
    char *strOut;
    int MAX_ELEMENT_LENGTH = 20;
    int MAX_ELEMENTS = 7;
    int length = 0;
    char *comment = "Output produced by medFlow2D";

    // Get file name parts
    char *baseFile = basename(strdup(fName));
    char *baseDir = dirname(strdup(fName));
    char archiveName[1024];
    char fullArchiveNAme[2048];

    // zip archive name
    sprintf(archiveName, "%s.zip", baseFile);
    sprintf(fullArchiveNAme, "%s/%s", baseDir, archiveName);

    //Guess size (over approximate)
    size_t length_approx = (size_t)(sim->ly*sim->lx*(MAX_ELEMENTS + MAX_ELEMENTS*MAX_ELEMENT_LENGTH));

    slog(5, SLOG_DEBUG, "Allocation output buffer of length %d...\n", length_approx);

    //allocate buffer
    strOut = calloc(length_approx, sizeof(char));

    slog(5, SLOG_DEBUG, "Filling up buffer (virtual %s)...\n", baseFile);

    // Write resolution
    length += bytes_added(sprintf(strOut+length, "%d\t%d\n", sim->lx, sim->ly));

    // Write output to buffer
    real ux, uy, fx, fy, rho, p, c;
    int iX, iY;
    for (iY=1; iY<=sim->ly; iY++) {
        for (iX=1; iX<=sim->lx; iX++) {
            if(sim->geometry[iX][iY]!=WALL && sim->geometry[iX][iY] != COAGWALL) {
                //computeMacros(sim->lattice[iX*step+1][iY*step+1].fPop, &rho, &ux, &uy);
                computeRho(sim->activeLattice[iX][iY].fPop, &rho);
                computeRho(sim->activePsLattice[iX][iY].fPop, &c);
                //p = ((rho-1.0)/3.0)*sim->dm/sim->dx/sim->dt/sim->dt;
                p = (rho - 1.0) / 3.0 * sim->dm / sim->dt / sim->dt;
                //p = (rho - 1.0) / 3.0 * sim->dx * sim->dx / sim->dt / sim->dt;
                //rho = rho *sim->dm/sim->dx/sim->dx/sim->dx;
                ux = sim->store[iX][iY].ux * sim->dx / sim->dt;
                uy = sim->store[iX][iY].uy * sim->dx / sim->dt;
                fx = sim->store[iX][iY].fx * sim->dm * sim->dx / sim->dt / sim->dt;
                fy = sim->store[iX][iY].fy * sim->dm * sim->dx / sim->dt / sim->dt;

                //uNorm = sqrt(ux*ux+uy*uy);
                length += bytes_added(
                        sprintf(strOut + length, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", ux, uy, rho, p, c - sim->rhoPSMean, fx,
                                fy));
                //fprintf(oFile, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", ux, uy, rho, sim->store[iX][iY].shear_max, sim->store[iX][iY].shear_ang, fx, fy);
            }
            else
                length += bytes_added(
                        sprintf(strOut + length, "%f\t%f\t%f\t%f\t%e\t%e\t%e\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        }
    }

    slog(5, SLOG_DEBUG, "Beginning compression of output data (%d bytes) to (%s)\n", length, fullArchiveNAme);

    // If file exists, remove it
    if(fileExists(fullArchiveNAme))
        remove(fullArchiveNAme);

    // Compression
    mz_bool status;
    status = mz_zip_add_mem_to_archive_file_in_place(fullArchiveNAme, baseFile, strOut, length + 1, comment, (short)strlen(comment), MZ_BEST_COMPRESSION);
    if (!status)
    {
        slog(0, SLOG_ERROR, "Failed to compress output file!\n");
        //return EXIT_FAILURE;
    }

    slog(5, SLOG_DEBUG, "Compression is done.\n");
    // Free up string buffer
    free(strOut);
}

int checkNanSparse(Simulation* sim, int step)
{
	int iX, iY;
	int maxX = floor(sim->lx/step);
	int maxY = floor(sim->ly/step);

	real rho;
	for (iY=0; iY<maxY; iY++)
		for (iX=0; iX<maxX; iX++) {
			computeRho(sim->activeLattice[iX*step+1][iY*step+1].fPop, &rho);
			if(isnan(rho))
				return -1;

			if(sim->coupled){
				computeRho(sim->activePsLattice[iX*step+1][iY*step+1].fPop, &rho);
				if(isnan(rho))
					return -1;
			}
		}

	return 0;
}

  // save the velocity field (norm) to disk
void saveVelPPM(Simulation *sim, char fName[], real uMax) {

    //printf("Creating image (%d x %d)\n", sim->lx, sim->ly);

    PPMImage *img = createPPM(sim->lx, sim->ly, uMax);

    int iX, iY;
    real ux, uy, uNorm, rho;
   
    for (iY=1; iY<=sim->ly; iY++) {
       for (iX=1; iX<=sim->lx; iX++) {
           if (sim->geometry[iX][iY] == WALL || sim->geometry[iX][iY] == COAGWALL || sim->geometry[iX][iY] == UNUSED)
               setPixelBlack(img, iX-1, iY-1);
           else {
               computeMacros(sim->activeLattice[iX][iY].fPop, &rho, &ux, &uy);
               uNorm = sqrt(ux*ux+uy*uy);
               setPixel(img, iX-1, iY-1, uNorm);
           }
       }
    }
    writePPM(fName, img);
    freePPM(img);
}

void saveVelPNG(Simulation* sim, char fName[], real uMax)
{
    int iXmax = sim->lx;
    int iYmax = sim->ly;

    char *pImage = (char *)malloc((size_t)(iXmax * 3 * iYmax));

    real ux, uy, uNorm, rho;

    int iX, iY;
    for(iY = 0; iY < iYmax; iY++)
    {
        for(iX = 0; iX < iXmax; iX++)
        {
            char *color = (char *)(pImage + (iX * 3) + (iY * iXmax * 3));

            if (sim->geometry[iX+1][iY+1] == WALL || sim->geometry[iX+1][iY+1] == COAGWALL || sim->geometry[iX+1][iY+1] == UNUSED) {
                color[0] = 0; color[1] = 0; color[2] = 0;
            }
            else {
                computeMacros(sim->activeLattice[iX+1][iY+1].fPop, &rho, &ux, &uy);
                uNorm = sqrt(ux*ux+uy*uy);

                real c = uNorm / uMax;

                if(c>1.0)
                    c = 1.0f;

                color[0] = (unsigned char)(255.0f*getJetR(c));
                color[1] = (unsigned char)(255.0f*getJetG(c));
                color[2] = (unsigned char)(255.0f*getJetB(c));

            }
        }
    }

    // Now save PNG
    size_t png_data_size = 0;
    void *pPNG_data = tdefl_write_image_to_png_file_in_memory_ex(pImage, iXmax, iYmax, 3, &png_data_size, 6, MZ_FALSE);
    if (!pPNG_data)
        slog(0, SLOG_ERROR, "Failed to save PNG image (%d)!\n", fName);
    else
    {
        FILE *pFile = fopen(fName, "wb");
        fwrite(pPNG_data, 1, png_data_size, pFile);
        fclose(pFile);
    }

    // mz_free() is by default just an alias to free() internally, but if you've overridden miniz's allocation funcs you'll probably need to call mz_free().
    mz_free(pPNG_data);
    free(pImage);
}

 // save the velocity field (norm) to disk
void saveRhoPPM(Simulation *sim, Node **lattice, char fName[], real RhoMean, real colorScale) {

    //printf("Creating image (%d x %d)\n", sim->lx, sim->ly);

    PPMImage *img = createPPM(sim->lx, sim->ly, 1.0);

    int iX, iY;
    real ux, uy, rho;
   
    for (iY=1; iY<=sim->ly; iY++) {
       for (iX=1; iX<=sim->lx; iX++) {
           if (sim->geometry[iX][iY] == WALL || sim->geometry[iX][iY] == COAGWALL || sim->geometry[iX][iY] == UNUSED)
               setPixelBlack(img, iX-1, iY-1);
           else {
               computeMacros(lattice[iX][iY].fPop, &rho, &ux, &uy);
               setPixel(img, iX - 1, iY - 1, ((rho - RhoMean) * colorScale) + RhoMean / 2.0);
           }
       }
    }
    writePPM(fName, img);
    freePPM(img);
}

void saveRhoPNG(Simulation *sim, Node **lattice, char fName[], real cMin, real cMax)
{
    int iXmax = sim->lx;
    int iYmax = sim->ly;

    char *pImage = (char *)malloc((size_t)(iXmax * 3 * iYmax));

    real rho;
    real scale = cMax - cMin;

    int iX, iY;
    for(iY = 0; iY < iYmax; iY++)
    {
        for(iX = 0; iX < iXmax; iX++)
        {
            char *color = (char *)(pImage + (iX * 3) + (iY * iXmax * 3));

            if (sim->geometry[iX+1][iY+1] == WALL || sim->geometry[iX+1][iY+1] == COAGWALL || sim->geometry[iX+1][iY+1] == UNUSED) {
                color[0] = 0; color[1] = 0; color[2] = 0;
            }
            else {
                computeRho(lattice[iX+1][iY+1].fPop, &rho);

                // Color scale
                real c = (rho-cMin)/scale;

                if(c>1.0)
                    c = 1.0f;
                if(c < 0.0)
                    c = 0.0f;

                color[0] = (unsigned char)(255.0f*getJetR(c));
                color[1] = (unsigned char)(255.0f*getJetG(c));
                color[2] = (unsigned char)(255.0f*getJetB(c));

            }
        }
    }

    // Now save PNG
    size_t png_data_size = 0;
    void *pPNG_data = tdefl_write_image_to_png_file_in_memory_ex(pImage, iXmax, iYmax, 3, &png_data_size, 6, MZ_FALSE);
    if (!pPNG_data)
        slog(0, SLOG_ERROR, "Failed to save PNG image (%d)!\n", fName);
    else
    {
        FILE *pFile = fopen(fName, "wb");
        fwrite(pPNG_data, 1, png_data_size, pFile);
        fclose(pFile);
    }

    // mz_free() is by default just an alias to free() internally, but if you've overridden miniz's allocation funcs you'll probably need to call mz_free().
    mz_free(pPNG_data);
    free(pImage);
}


 // save the velocity field (norm) to disk
void saveGeomImg(Simulation* sim,char fName[], int maxNodeType) {

    PPMImage *img = createPPM(sim->lx, sim->ly, (float)maxNodeType);

    int iX, iY;
   
    for (iY=1; iY<=sim->ly; iY++) 
       for (iX=1; iX<=sim->lx; iX++) 
           setPixel(img, iX-1, iY-1, sim->geometry[iX][iY]);
    
    
    writePPM(fName, img);
    freePPM(img);
}


 // save the velocity field (norm) to disk
void saveRho(Simulation* sim, Node** lattice, char fName[]) {
    FILE* oFile = fopen(fName, "w");
    int iX, iY;
    real ux, uy, rho;
    for (iY=1; iY<=sim->ly; iY++) {
       for (iX=1; iX<=sim->lx; iX++) {
   		   computeMacros(lattice[iX][iY].fPop, &rho, &ux, &uy);
   		   fprintf(oFile, "%f ", rho);
       }
       fprintf(oFile, "\n");
    }
    fclose(oFile);
}
  // save one lattice population to disk
void saveF(Simulation* sim, Node** lattice, int iPop, char fName[]) {
    FILE* oFile = fopen(fName, "w");
    int iX, iY;

    for (iY=1; iY<=sim->ly; iY++) {
       for (iX=1; iX<=sim->lx; iX++) {
           real f = lattice[iX][iY].fPop[iPop];
           fprintf(oFile, "%f ", f);
       }
       fprintf(oFile, "\n");
    }
    fclose(oFile);
}


void saveScalar(char fName[], int lx, int ly, real **field)
{
    FILE* oFile = fopen(fName, "w");

    int iX, iY;

    for (iY=1; iY<=ly; iY++) {
       for (iX=1; iX<=lx; iX++) {           
           fprintf(oFile, "%f ", field[iX][iY]);
       }
       fprintf(oFile, "\n");
    }
    fclose(oFile);
}

  // save the velocity field (norm) to disk
void saveScalarImg(char fName[], int lx, int ly, real uMax, real** array) {

    PPMImage *img = createPPM(lx, ly, uMax);

    int iX, iY;
   
    for (iY=1; iY<=ly; iY++)
       for (iX=1; iX<=lx; iX++)
    	   setPixel(img, iX-1, iY-1, array[iX][iY]);

    writePPM(fName, img);
    freePPM(img);
}
/*
PPM read/write code from:
http://stackoverflow.com/questions/2693631/read-ppm-file-and-store-it-in-an-array-coded-with-c
*/

static inline real clamp(real value){
  if(value>1.0)
    return 1.0;
  else if(value<0.0)
    return 0.0;
  else
    return value;
}

static inline real getJetR(real gray){
  real fGray = 4.0*gray;
  return clamp(min(fGray - 1.5, -fGray + 4.5));
}

static inline real getJetG(real gray){
  real fGray = 4.0*gray;
  return clamp(min(fGray - 0.5, -fGray + 3.5));
}

static inline real getJetB(real gray){
  real fGray = 4.0*gray;
  return clamp(min(fGray + 0.5, -fGray + 2.5));
}

PPMImage *createPPM(const int w, const int h, const real maxInt) 
{
	PPMImage *img;
	img = (PPMImage *)malloc(sizeof(PPMImage));

    if (!img) {
        slog(1, SLOG_ERROR , "Insufficient memory for (%d x %d) image.\n", w, h);
        exit(1);
    }

    img->x = w; img->y = h; img->maxIntensity = maxInt;
	img->data = (PPMPixel*)calloc(w * h ,sizeof(PPMPixel));

    return img;
}

void setPixel(PPMImage* img, const int x, const int y, const real color) {
	// Convert intensity to [0,1] scale
	real c = color / img->maxIntensity;
	if(c>1.0){
		//printf("Warning,color intensity overflow! (%f of %f)\n", color, img->maxIntensity);
		c = 1.0f;
	}
	
	int loc = y*img->x+x;
	img->data[loc].r = (unsigned char)(255.0f*getJetR(c));
	img->data[loc].g = (unsigned char)(255.0f*getJetG(c));
	img->data[loc].b = (unsigned char)(255.0f*getJetB(c));
}

void setPixelBlack(PPMImage* img, const int x, const int y) {
    int loc = y*img->x+x;
    img->data[loc].r = (unsigned char)0;
    img->data[loc].g = (unsigned char)0;
    img->data[loc].b = (unsigned char)0;
}

PPMImage *readPPM(const char *filename) 
{
    char buff[16];
    PPMImage *img;
    FILE *fp;
    int c, rgb_comp_color;
  
    fp = fopen(filename, "rb");
    if (!fp) {
        slog(0, SLOG_ERROR, "Error opening file: '%s'\n", filename);
        exit(1);
    }

    if (!fgets(buff, sizeof(buff), fp)) {
        perror(filename);
        exit(1);
    }

    if (buff[0] != 'P' || buff[1] != '6') {
        slog(0, SLOG_ERROR, "Wrong file format! (need 'P6')\n");
        exit(1);
    }

    img = (PPMImage *)malloc(sizeof(PPMImage));
    if (!img) {
        slog(0, SLOG_ERROR, "Insufficient memory!\n");
        exit(1);
    }

    c = getc(fp);
    while (c == '#') {
        while (getc(fp) != '\n') ;
        c = getc(fp);
    }

    ungetc(c, fp);
    
    if (fscanf(fp, "%d %d", &img->x, &img->y) != 2){
        slog(0, SLOG_ERROR, "Image size error ('%s')\n", filename);
        exit(1);
    }

    if (fscanf(fp, "%d", &rgb_comp_color) != 1)  {
        slog(0, SLOG_ERROR, "Wrong RGB value ('%s')\n", filename);
        exit(1);
    }

    if (rgb_comp_color!= CH_BITS) {
        slog(0, SLOG_ERROR, "'%s' is not an8-bit image!\n", filename);
        exit(1);
    }

    while (fgetc(fp) != '\n') ;
	
	img->data = (PPMPixel*)malloc(img->x * img->y * sizeof(PPMPixel));

    if (!img) {
        fprintf(stderr, "Insufficient memory.\n");
        exit(1);
    }

    if (fread(img->data, 3 * img->x, img->y, fp) != img->y)  {
        slog(0, SLOG_ERROR, "Reading error: '%s'\n", filename);
        exit(1);
    }

    fclose(fp);
    return img;
}

void writePPM(const char *filename, PPMImage *img)
{
    FILE *fp;

    fp = fopen(filename, "wb");
    if (!fp) {
        slog(0, SLOG_ERROR, "Unable to open file '%s'\n", filename);
        exit(1);
    }

    fprintf(fp, "P6\n");

    fprintf(fp, "# Created by medFlow2D. Copyright 2015, Gabor Zavodszky.\n");

    fprintf(fp, "%d %d\n",img->x,img->y);

    fprintf(fp, "%d\n",CH_BITS);

    fwrite(img->data, 3 * img->x, img->y, fp);
    fclose(fp);
}

void freePPM(PPMImage *img)
{
    free(img->data);
    free(img);
}

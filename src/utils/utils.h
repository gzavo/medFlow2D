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


#ifndef UTILS_H
#define UTILS_H

#include "../common.h"
#include "../lbm/lb.h"

static const char* cmdFile = ".command";

int checkCommand(Simulation *sim, const char * fileName);

/*
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
*/

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

inline real interpolate1(real x1, real x2, real xi, real y1, real y2);
real poiseuilleProfile(int r, real uMax, int R);

typedef struct
{
	int numVal;
	real *scales;
	real *times;
	int time_current;	// The currently active discrete time
	real time_max;
	int currentCycle;
} transientBC;

real inline getScale(transientBC *bc, real phys_time);
transientBC *loadTransientBC(const char *fileName);

void printBC(transientBC *bc);
void printSimDetails(Simulation* sim);
void fprintSimDetails(Simulation* sim, char *fileName);

// ProgressBar
void loadBar(int x, int n, int r, int w, int cells);

// Data saving routines

void saveVel(Simulation* sim, char fName[]);
void saveRho(Simulation* sim, Node **lattice, char fName[]);
void saveVelPPM(Simulation *sim, char fName[], real uMax);
void saveVelPNG(Simulation* sim, char fName[], real uMax);
void saveRhoPPM(Simulation *sim, Node **lattice, char fName[], real RhoMean, real colorScale);
void saveRhoPNG(Simulation *sim, Node **lattice, char fName[], real cMin, real cMax);
void saveFullSparse(Simulation* sim, char fName[], int step);
void saveFullSparseZ(Simulation* sim, char fName[], int step);
void saveFull(Simulation* sim, char fName[]);
void saveFullZ(Simulation* sim, char fName[]);
void saveF(Simulation* sim, Node** lattice, int iPop, char fName[]);
void saveScalar(char fName[], int lx, int ly, real **field);
void saveScalarImg(char fName[], int lx, int ly, real uMax, real** array);
void saveGeomImg(Simulation* sim,char fName[], int maxNodeType) ;

int checkNanSparse(Simulation* sim, int step);

typedef struct
{
	int x, y;
} coord;

typedef enum {
	NOCMD=0,
	SAVEA,
	SAVEF,
	SAVEV,
	SAVED,
	SAVEP
} Signals;

typedef struct {
	char *commandStr;
	Signals sign;
} Command;

char *baseName(char const *path);

// PPM stuff

typedef struct
{
    unsigned char r,g,b;
} PPMPixel;

typedef struct
{
    int x, y;
    float maxIntensity;		// Maximal intensity value
    PPMPixel *data;
} PPMImage;

// 1 byte/colorchannel
#define CH_BITS 255

// File handling routines
int dirExists(const char *path);
int fileExists(const char *fname);
int createDir(const char *path);
int copyFile(char *old_filename, char  *new_filename);

void saveCheckPoint(Simulation* sim, const char *fName);
void loadCheckPoint(Simulation* sim, const char *fName);

int checkPSSanity(Simulation *sim);

PPMImage *createPPM(const int w, const int h, const real maxInt);
void setPixel(PPMImage* img, const int x, const int y, const real color);
void setPixelBlack(PPMImage* img, const int x, const int y);
PPMImage *readPPM(const char *filename);
void writePPM(const char *filename, PPMImage *img);
void freePPM(PPMImage *img);

#endif

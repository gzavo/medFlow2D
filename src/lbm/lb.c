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

#define THREAD_PARALLEL     // switch between thread and OpenMP based parallelization

#include <stdlib.h>
#include <math.h>
#include "lb.h"
#include "../utils/tinycthread.h"

//#include <stdio.h>
//#include <float.h>
//#include <pthread.h>
//#include "assert.h"
//#include "../utils/utils.h"

// For C99 inline handling
extern real inline computeEquilibrium(int iPop, real rho, real ux, real uy, real uSqr);

extern void inline computeMRTEquilibrium(real *m_eq, real rho, real ux, real uy);

extern void inline computeMacros(real *f, real *rho, real *ux, real *uy);

extern real inline computeForce(int iPop, real rho, real omega, real ux, real uy, real Fx, real Fy, real Fu);

extern void inline porousForcing(real *fx, real *fy, const real c0, const real c1, const real ux, const real uy);

extern void inline mrt2bgk(real *mPop, real *fPop);

extern void inline mrtCollision(real *m, real *mEq, void *selfData);

extern void inline bgk2mrt(real *fPop, real *mPop);

extern void inline computeRho(real *f, real *rho);

extern void inline calculateShearStressMax(NodeStore *store, real *shearMax, real *shearAngle);

extern real inline computePSEquilibrium(int iPop, real rho, real ux, real uy);



/* Global state modifiers */
/**************************/
//TODO They block parallel execution? - most likely
/*bool computeStressTensor = false;
real margForceRatio = 1.0;

void setComputeStressTensor(bool value)
{
   computeStressTensor = value;
}

inline bool getComputeStressTensor()
{ return computeStressTensor; }

real getMargForceRatio ()
{ return margForceRatio; }

void setMargForceRatio(real mFR)
{ margForceRatio = mFR; }

*/
//static const int downAA[4] = {1,2,3,4};
//static const int upAA[4] = {5,6,7,8};

/* MRT arrays arrays */
/**************/
const real M[9][9] = {
		{ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.},
		{-4., -1.,  2., -1.,  2., -1.,  2., -1.,  2.},
		{ 4., -2.,  1., -2.,  1., -2.,  1., -2.,  1.},
		{ 0.,  1.,  1.,  0., -1., -1., -1.,  0.,  1.},
		{ 0., -2.,  1.,  0., -1.,  2., -1.,  0.,  1.},
		{ 0.,  0., -1., -1., -1.,  0.,  1.,  1.,  1.},
		{ 0.,  0., -1.,  2., -1.,  0.,  1., -2.,  1.},
		{ 0.,  1.,  0., -1.,  0.,  1.,  0., -1.,  0.},
		{ 0.,  0., -1.,  0.,  1.,  0., -1.,  0.,  1.}
};

const real Minv[9][9] = {
		{ 1./9, -1./ 9,  1./ 9,  0   ,  0    ,  0   ,  0    ,  0   ,  0   },
		{ 1./9, -1./36, -1./18,  1./6, -1./ 6,  0   ,  0    ,  1./4,  0   },
		{ 1./9,  1./18,  1./36,  1./6,  1./12, -1./6, -1./12,  0   , -1./4},
		{ 1./9, -1./36, -1./18,  0   ,  0    , -1./6,  1./ 6, -1./4,  0   },
		{ 1./9,  1./18,  1./36, -1./6, -1./12, -1./6, -1./12,  0   ,  1./4},
		{ 1./9, -1./36, -1./18, -1./6,  1./6 ,  0   ,  0    ,  1./4,  0   },
		{ 1./9,  1./18,  1./36, -1./6, -1./12,  1./6,  1./12,  0   , -1./4},
		{ 1./9, -1./36, -1./18,  0   , 0     ,  1./6, -1./ 6, -1./4,  0   },
		{ 1./9,  1./18,  1./36,  1./6,  1./12,  1./6,  1./12,  0   ,  1./4}
};
/*
static const int M[9][9] = {
   {1, 1, 1, 1, 1, 1, 1, 1, 1},
  {-4,-1,-1,-1,-1, 2, 2, 2, 2}, 
   {4,-2,-2,-2,-2, 1, 1, 1, 1},
   {0, 1, 0,-1, 0, 1,-1,-1, 1},
   {0,-2, 0, 2, 0, 1,-1,-1, 1},
   {0, 0, 1, 0,-1, 1, 1,-1,-1},
   {0, 0,-2, 0, 2, 1, 1,-1,-1},
   {0, 1,-1, 1,-1, 0, 0, 0, 0},
   {0, 0, 0, 0, 0, 1,-1, 1,-1}
};
*/
/*
static const int Minv[9][9] = {
    {1./9, -1./9,   1./9,      0,      0,     0,      0,     0,     0}, 
    {1./9, -1./36, -1./18,  1./6,  -1./6,     0,      0,  1./4,     0}, 
    {1./9, -1./36, -1./18,     0,      0,  1./6,  -1./6, -1./4,     0},  
    {1./9, -1./36, -1./18, -1./6,   1./6,     0,      0,  1./4,     0}, 
    {1./9, -1./36, -1./18,     0,      0, -1./6,   1./6, -1./4,     0}, 
    {1./9,  1./18,  1./36,  1./6,  1./12,  1./6,  1./12,     0,  1./4}, 
    {1./9,  1./18,  1./36, -1./6, -1./12,  1./6,  1./12,     0, -1./4}, 
    {1./9,  1./18,  1./36, -1./6, -1./12, -1./6, -1./12,     0,  1./4}, 
    {1./9,  1./18,  1./36,  1./6,  1./12, -1./6, -1./12,     0, -1./4} 
};
*/

//								 *	    *	       * 	      *  <- free to set these
//						  rho    e      eps    jx  qx     jy  qy   pxx  pxy
const real S[9] = {0  ,   1.001, 1.001, 0 , 1.001, 0 , 1.001, 0  , 0  };
//static const real S[9] = {0  ,   1.0, 1.0, 0 , 1.0, 0 , 1.0, 0  , 0  };




Dynamics *createDynamics() {

	return NULL;
}


/* struct Node, methods                                          */
/*****************************************************************/

// initialize a node to default value
void constructNode(Node *node) {
	int iPop;

	for (iPop = 0; iPop < 9; ++iPop) {
		node->fPop[iPop] = 0.;
	}
	node->dynamics = 0;
}

// initialize a node to its local equilibrium term
void iniEquilibrium(Simulation *sim, int x, int y, real rho, real ux, real uy) {
	int iPop;
	real uSqr = ux * ux + uy * uy;

	for (iPop = 0; iPop < 9; ++iPop) {
		sim->latticeA[x][y].fPop[iPop] = computeEquilibrium(iPop, rho, ux, uy, uSqr);
	}
	if(sim->doubleBuffered)
		cpyDistribution(sim, sim->latticeA, sim->latticeB);
}

// initialize a node to its local equilibrium term
void iniPSEquilibrium(Simulation *sim, int x, int y, real rho, real ux, real uy) {
	int iPop;
	real uSqr = ux * ux + uy * uy;

	for (iPop = 0; iPop < 9; ++iPop) {
		sim->psLatticeA[x][y].fPop[iPop] = computeEquilibrium(iPop, rho, ux, uy, uSqr);
	}
	if(sim->doubleBuffered)
		cpyDistribution(sim, sim->psLatticeA, sim->psLatticeB);
}

void flipActiveLattice(Simulation *sim)
{
	if(sim->doubleBuffered)
	{
		if(sim->activeLattice == sim->latticeA) {
			sim->activeLattice = sim->latticeB;
			sim->activePsLattice = sim->psLatticeB;
		}
		else {
			sim->activeLattice = sim->latticeA;
			sim->activePsLattice = sim->psLatticeA;
		}
	}
}

/* struct Simulation, methods                                    */
/*****************************************************************/

// allocate memory for a full simulation ("constructor")
void constructSim(Simulation *sim, int lx, int ly, int doDoubleBuffer, int coupled, int thrombus) {
	sim->lx = lx;
	sim->ly = ly;

	sim->doubleBuffered = doDoubleBuffer;
	sim->coupled = coupled;
	sim->thrombus = thrombus;

	// Allocate memory for Buffer A and geometry
	sim->memoryChunkA = (Node *) calloc((lx + 2) * (ly + 2), sizeof(Node));
	sim->latticeA = (Node **) calloc(lx + 2, sizeof(Node *));

	// Set active lattice to A
	sim->activeLattice = sim->latticeA;

	sim->geometryMemoryChunk = (NodeType *) calloc((lx + 2) * (ly + 2), sizeof(NodeType));
	sim->geometry = (NodeType **) calloc(lx + 2, sizeof(NodeType *));

	// Allocate tmpLattice anyhow
	sim->storeMemoryChunk = (NodeStore *) calloc((lx + 2) * (ly + 2), sizeof(NodeStore));
	sim->store = (NodeStore **) calloc(lx + 2, sizeof(NodeStore *));

	sim->psMemoryChunkA = (Node *) calloc((lx + 2) * (ly + 2), sizeof(Node));
	sim->psLatticeA = (Node **) calloc(lx + 2, sizeof(Node *));

	// Set active tmpLattice
	sim->activePsLattice = sim->psLatticeA;


	int iX, iY;
	for (iX = 0; iX < lx + 2; ++iX) {

		sim->latticeA[iX] = sim->memoryChunkA + iX * (ly + 2);
		sim->geometry[iX] = sim->geometryMemoryChunk + iX * (ly + 2);

		// Allocate anyhow
		//if(coupled>0) {
		sim->psLatticeA[iX] = sim->psMemoryChunkA + iX * (ly + 2);
		sim->store[iX] = sim->storeMemoryChunk + iX * (ly + 2);
		//}

		for (iY = 0; iY < ly + 2; ++iY) {
			constructNode(&(sim->latticeA[iX][iY]));
			setGeometry(sim, iX, iY, FLUID);
			//if(coupled>0)
			constructNode(&(sim->psLatticeA[iX][iY]));
		}
	}

	if (sim->doubleBuffered) {
		// Allocate memory for Buffer B
		sim->memoryChunkB = (Node *) calloc((lx + 2) * (ly + 2), sizeof(Node));
		sim->latticeB = (Node **) calloc(lx + 2, sizeof(Node *));

		sim->psMemoryChunkB = (Node *) calloc((lx + 2) * (ly + 2), sizeof(Node));
		sim->psLatticeB = (Node **) calloc(lx + 2, sizeof(Node *));

		int iX, iY;
		for (iX = 0; iX < lx + 2; ++iX) {

			sim->latticeB[iX] = sim->memoryChunkB + iX * (ly + 2);

			sim->psLatticeB[iX] = sim->psMemoryChunkB + iX * (ly + 2);

			for (iY = 0; iY < ly + 2; ++iY) {
				constructNode(&(sim->latticeB[iX][iY]));
				constructNode(&(sim->psLatticeB[iX][iY]));
			}
		}

	}

}

// free the memory for the simulation ("destructor") - does not free simulation dynamics
void destructSim(Simulation *sim) {
	free(sim->latticeA);
	free(sim->memoryChunkA);
	free(sim->store);
	free(sim->storeMemoryChunk);

	//Always allocated as coupled
	free(sim->psLatticeA);
	free(sim->psMemoryChunkA);

	if (sim->doubleBuffered) {
		free(sim->latticeB);
		free(sim->memoryChunkB);
		free(sim->psLatticeB);
		free(sim->psMemoryChunkB);
	}
}

// copy values from main lattice
void cpyDistribution(Simulation *sim, Node **latticeA, Node **latticeB) {
	int i, j;
	for (i = 0; i <= sim->lx; i++)
		for (j = 0; j <= sim->ly; j++) {
			int k;
			for (k = 0; k < 9; k++) {
				latticeB[i][j].fPop[k] = latticeA[i][j].fPop[k];
			}
		}
}


// specify the dynamics for a given lattice site
void setDynamics(Simulation *sim, int x, int y, Dynamics *dyn) {
	sim->latticeA[x][y].dynamics = dyn;
	if(sim->doubleBuffered)
		sim->latticeB[x][y].dynamics = dyn;
}

void setPSDynamics(Simulation *sim, int x, int y, Dynamics *dyn) {
	sim->psLatticeA[x][y].dynamics = dyn;
	if(sim->doubleBuffered)
		sim->psLatticeB[x][y].dynamics = dyn;
}


Dynamics *getDynamics(Node *cell) {
	return cell->dynamics;
}

void setGeometry(Simulation *sim, int iX, int iY, NodeType type) {
	sim->geometry[iX][iY] = type;
}

NodeType getGeometry(Simulation *sim, int iX, int iY) {
	return sim->geometry[iX][iY];
}


// apply collision step to a lattice node (and simulate
//   virtual dispatch)
inline static void collideNodeAA(Node *node, NodeStore *store) {
	node->dynamics->dynamicsFun(node->fPop, node->dynamics->selfData, store);
}

// Parallelise with pthreads
int collideWorker(void *arg) {

	ParallelData *p = (ParallelData *) arg;
	Simulation *sim = p->sim;
	Node **lattice = p->nodes;

	int iX, iY;
	for (iX = p->fromIdx; iX <= p->toIdx; iX++) {
		for (iY = 1; iY <= sim->ly; iY++) {
            if (getGeometry(sim, iX, iY) != UNUSED)
			    collideNodeAA(&(lattice[iX][iY]), &sim->store[iX][iY]);
		}
	}

	return 0;
}

// apply collision step to all lattice nodes
inline void collideAA(Simulation *sim, Node **lattice, int numWorkers) {

#ifdef THREAD_PARALLEL

    thrd_t *t = (thrd_t*)malloc(numWorkers * sizeof(thrd_t));
    ParallelData *ptData = (ParallelData *)malloc(numWorkers * sizeof(ParallelData));

    int numPerProcess = (int)floor((real) sim->lx / numWorkers);

    int i;
    for(i = 0; i < numWorkers; i++) {

        ptData[i].sim = sim;
        ptData[i].nodes = lattice;
        ptData[i].id = i;

        ptData[i].fromIdx = i * numPerProcess + 1; //1 is the first cell (see ghost cells)

        if (i == numWorkers - 1)  //last thread does some overwork
            ptData[i].toIdx = sim->lx;
        else
            ptData[i].toIdx = (i + 1) * numPerProcess;

        thrd_create(&t[i], &collideWorker, (void*)&ptData[i]);
    }

    for(i=0; i < numWorkers; i++) {
        thrd_join(t[i], NULL);
    }

    free(t);
    free(ptData);
#else

	int iX, iY;
	int maxX = sim->lx;
	int maxY = sim->ly;
	NodeStore **ns = sim->store;
	//#pragma omp parallel for private(iX, iY) // shared(lattice, sim)

#pragma omp parallel for shared(sim, lattice, ns, maxX, maxY) private(iX, iY)
		for (iX = 1; iX <= maxX; iX++) {
			for (iY = 1; iY <= maxY; iY++) {
				if (getGeometry(sim, iX, iY) != UNUSED)
					collideNodeAA(&(lattice[iX][iY]), &ns[iX][iY]);
			}
		}
#endif

}

void calculateShearStressMaxDomain(Simulation *sim, real **shearMax, real **shearAngle) {
	int iX, iY;

	for (iX = 1; iX <= sim->lx; iX++) {
		for (iY = 1; iY <= sim->ly; iY++)
			calculateShearStressMax(&sim->store[iX][iY], &shearMax[iX][iY], &shearAngle[iX][iY]);
	}
}

inline void calculateShearStressMax(NodeStore *store, real *shearMax, real *shearAngle) {
	real sigma[3];
	real mu2 = 2.0 * store->mu;
	sigma[0] = mu2 * store->Sij[0];
	sigma[1] = mu2 * store->Sij[1];
	sigma[2] = mu2 * store->Sij[2];

	// Maximal shear stress and direction from Mohr's diagram
	real sqrDist = (sigma[0] - sigma[2]) / 2.0;
	real r = sqrt(sqrDist * sqrDist + sigma[1] * sigma[1]);

	//real s_A = (s[0]+s[2])/2;
	//real ang_P = atan(-s[1]/(s[1]-s_A))/2.0;  // Primal stress direction

	real ang_P = atan2(sigma[1], sqrDist) / 2.0;
	//real ang_P = atan(sigma[1]/sqrDist);

	// TODO: check! i'm not sure if its correct!!!
	real ang_S = M_PI / 2.0 - ang_P; //+ copysign(M_PI/2.0, sigma[1]);// - M_PI/2.0; //+ copysign(M_PI/4.0, sigma[1]);              // 45deg - s_P

	// For y aligned pipe:
	//real ang_S = ang_P - copysign(M_PI/2.0, sigma[0]*sigma[2]);

	//store->shear_max = r;
	//store->shear_ang = ang_S;

	*shearMax = r;
	*shearAngle = ang_S;
}


inline void propagateAA(Simulation *sim, Node **lattice) {
	int iX, iY;
	int lx = sim->lx;
	int ly = sim->ly;

	//Propagate downwards (towards positive x and positive y)
	for (iY = 1; iY <= ly; iY++) {
		for (iX = 1; iX <= lx; iX++) {
			if (getGeometry(sim, iX, iY) != UNUSED){
				lattice[iX + c[2][0]][iY + c[2][1]].fPop[2] = lattice[iX][iY].fPop[2];
				lattice[iX + c[3][0]][iY + c[3][1]].fPop[3] = lattice[iX][iY].fPop[3];
				lattice[iX + c[4][0]][iY + c[4][1]].fPop[4] = lattice[iX][iY].fPop[4];
				lattice[iX + c[5][0]][iY + c[5][1]].fPop[5] = lattice[iX][iY].fPop[5];
			}
		}
	}

	//Propagate upwards (towards negative x and negative y)
	for (iY = ly; iY >= 1; iY--) {
		for (iX = lx; iX >= 1; iX--) {
			if (getGeometry(sim, iX, iY) != UNUSED){
				lattice[iX + c[6][0]][iY + c[6][1]].fPop[6] = lattice[iX][iY].fPop[6];
				lattice[iX + c[7][0]][iY + c[7][1]].fPop[7] = lattice[iX][iY].fPop[7];
				lattice[iX + c[8][0]][iY + c[8][1]].fPop[8] = lattice[iX][iY].fPop[8];
				lattice[iX + c[1][0]][iY + c[1][1]].fPop[1] = lattice[iX][iY].fPop[1];
			}
		}
	}



}

// implement periodic boundary conditions (to be called after
//   the propagation step)
void makePeriodic(Simulation *sim, Node **lattice) {
	int lx = sim->lx;
	int ly = sim->ly;
	Node **lat = lattice;

	int iX, iY;
	for (iX = 1; iX <= lx; ++iX) {
		lat[iX][ly].fPop[4] = lat[iX][0].fPop[4];
		lat[iX][ly].fPop[7] = lat[iX][0].fPop[7];
		lat[iX][ly].fPop[8] = lat[iX][0].fPop[8];

		lat[iX][1].fPop[2] = lat[iX][ly + 1].fPop[2];
		lat[iX][1].fPop[5] = lat[iX][ly + 1].fPop[5];
		lat[iX][1].fPop[6] = lat[iX][ly + 1].fPop[6];
	}

	for (iY = 1; iY <= ly; ++iY) {
		lat[1][iY].fPop[1] = lat[lx + 1][iY].fPop[1];
		lat[1][iY].fPop[5] = lat[lx + 1][iY].fPop[5];
		lat[1][iY].fPop[8] = lat[lx + 1][iY].fPop[8];

		lat[lx][iY].fPop[3] = lat[0][iY].fPop[3];
		lat[lx][iY].fPop[6] = lat[0][iY].fPop[6];
		lat[lx][iY].fPop[7] = lat[0][iY].fPop[7];
	}

	lat[1][1].fPop[5] = lat[lx + 1][ly + 1].fPop[5];
	lat[lx][1].fPop[6] = lat[0][ly + 1].fPop[6];
	lat[lx][ly].fPop[7] = lat[0][0].fPop[7];
	lat[1][ly].fPop[8] = lat[lx + 1][0].fPop[8];
}



/* some free helper functions                                    */
/*****************************************************************/

// compute density and velocity from the f's
void inline computeMacros(real *f, real *rho, real *ux, real *uy) {
	real upperLine = f[4] + f[3] + f[2];
	real mediumLine = f[5] + f[0] + f[1];
	real lowerLine = f[6] + f[7] + f[8];
	*rho = upperLine + mediumLine + lowerLine;
	*ux = (f[1] + f[2] + f[8] - (f[4] + f[5] + f[6])) / (*rho);
	*uy = (lowerLine - upperLine) / (*rho);
}

// compute rho only for passive scalar arrays (ux,uy comes from main lattice)
void inline computeRho(real *f, real *rho) {
	*rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
}

void equBC(real *fPop, void *selfData, NodeStore *store) {
	int iPop;
	real rho = ((real *) selfData)[0];
	real ux = ((real *) selfData)[1];
	real uy = ((real *) selfData)[2];

	store->ux = ux;
	store->uy = uy;

	real uSqr = ux * ux + uy * uy;

	for (iPop = 0; iPop < 9; ++iPop) {
		fPop[iPop] = computeEquilibrium(iPop, rho, ux, uy, uSqr);
	}
}

/* BGK Dynamics for porous media                                   */
/*******************************************************************/

// compute local equilibrium from rho and u
real inline computeForce(const int iPop, const real rho, const real omega, const real ux, const real uy, const real Fx, const real Fy, const real Fu) {
	real c_F = c[iPop][0] * Fx + c[iPop][1] * Fy;
	real c_u = c[iPop][0] * ux + c[iPop][1] * uy;
	return rho * t[iPop] * (1 - omega * 0.5) * (3. * c_F + 4.5 * c_u * c_F - 1.5 * Fu);
}

void inline porousForcing(real *fx, real *fy, const real c0, const real c1, const real ux, const real uy)
{
    *fx = -(c0 * ux + c1 * ux * fabs(ux));
    *fy = -(c0 * uy + c1 * uy * fabs(uy));
}

// bgk collision term
void inline bgkForced(real *fPop, void *selfData, NodeStore *store) {

	porousData *pd = (porousData *) selfData;

	real omega = pd->fp.omega;
	//real c0 = pd->c0;        // Darcy linear term
	//real c1 = pd->c1;        // Darcy quadratic term

	real rho, ux, uy, uSqr, Fx, Fy;

    real forcePop[9];
    real f_eq[9];

    computeMacros(fPop, &rho, &ux, &uy);

	// Compute modified macroscopic velocity

	//real u_mag = sqrt(ux * ux + uy * uy);
	//real c = c0+sqrt(c0*c0 + c1*sqrt(uSqr));
	// The negative sign is for the resistance opposes the velocity

	//Fx = -(c0*ux + c1*ux*u_mag);
	//Fy = -(c0*uy + c1*uy*u_mag);

    // Calculate forcing term
    porousForcing(&Fx, &Fy, pd->c0, pd->c1, ux, uy);

	real Fu = Fx * ux + Fy * uy;

	int iPop;
	for (iPop = 0; iPop < 9; iPop++)
		forcePop[iPop] = computeForce(iPop, rho, omega, ux, uy, Fx, Fy, Fu);

	// Modify macroscopic velocity
	ux += Fx / 2.0 / rho;
	uy += Fy / 2.0 / rho;

	// Store macroscopic values
	store->ux = ux;
	store->uy = uy;

	store->fx = Fx;
	store->fy = Fy;

	if (pd->fp.computeStressTensor) {
		store->Sij[0] = 0;
		store->Sij[1] = 0;
		store->Sij[2] = 0;
	}

	uSqr = ux * ux + uy * uy;

	// Compute equilibrium and non-equilibrium part of the stress tensor
	for (iPop = 0; iPop < 9; iPop++) {
		f_eq[iPop] = computeEquilibrium(iPop, rho, ux, uy, uSqr);

		if (pd->fp.computeStressTensor) {
			real f_neq = fPop[iPop] - f_eq[iPop];
			store->Sij[0] += f_neq * c[iPop][0] * c[iPop][0];
			store->Sij[1] += f_neq * c[iPop][0] * c[iPop][1];
			store->Sij[2] += f_neq * c[iPop][1] * c[iPop][1];
		}
	}

	// Compute the strain rate tensor from the stress tensor
	if (pd->fp.computeStressTensor) {
		//real sigmaFact =  - omega / (1./(3.0*omega) + 0.1666666);    //TODO: i'm not sure about this...
		real sigmaFact = -omega / (2. / 3. * rho);
		store->Sij[0] *= sigmaFact;
		store->Sij[1] *= sigmaFact;
		store->Sij[2] *= sigmaFact;
		store->mu = rho * (1.0 / (3.0 * omega) - 1.0 / 6.0);
	}

	// Recompute local omega for non-Newtonian dynamics
	if (pd->fp.useNonNewt) {
		NonNewtParam *np = pd->fp.np;

		real shearRate = sqrt(store->Sij[0] * store->Sij[0] +
				2 * store->Sij[1] * store->Sij[1] +
				store->Sij[2] * store->Sij[2]);

		real nu = np->nu_inf + (np->nu_0 - np->nu_inf) * pow(1 + pow(np->lambda * shearRate, np->a), (np->n - 1) / np->a);
		omega = 1.0 / (3.0 * nu + 0.5);
        pd->fp.omega = omega;
	}

    // Apply Smagorinsky turbulence model
    if(pd->fp.useSmagorinsky) {
        real tau0 = 1./omega;
        real sqrtStress = sqrt(store->Sij[0]*store->Sij[0]+2.0*store->Sij[1]*store->Sij[1]+store->Sij[2]*store->Sij[2]);
        real tau_t = 0.5*(sqrt(tau0*tau0 + 18.0*pd->fp.smagorinskyC*pd->fp.smagorinskyC*sqrtStress)-tau0);

        omega = 1.0/(tau0+tau_t);
    }

    real omega_min = 1.0 - omega;

    // Compute forced LBM collision result
	for(iPop=0; iPop < 9; iPop++)
	{
		fPop[iPop] *= omega_min;
		fPop[iPop] += omega * f_eq[iPop] + forcePop[iPop];
	}
}

/* Passive scalar dynamics										   */
/*******************************************************************/

real inline computePSEquilibrium(int iPop, real rho, real ux, real uy) {
	return rho * t[iPop] * (1. + 3. * (c[iPop][0] * ux + c[iPop][1] * uy));
}

// Passive scalar BGK
void inline bgkPS(real *fPop, void *selfData, NodeStore *store) {
	plateletData *ps = (plateletData *) selfData;
	real omega = ps->omega;
	real rho;

	computeRho(fPop, &rho);

	int iPop;
	for (iPop = 0; iPop < 9; iPop++) {
		fPop[iPop] *= (1 - omega);
		fPop[iPop] += omega * computePSEquilibrium(iPop, rho, store->ux, store->uy);
	}

}

void inline bgkForcedPS(real *fPop, void *selfData, NodeStore *store) {
	plateletData *ps = (plateletData *) selfData;
	real omega = ps->omega;
	real rho;

	computeRho(fPop, &rho);


	// Computation of shear stress maxes based on Mohr's circles, and on Kronenburg 2013 arxiv preprint
	real a, b, c;
	a = store->Sij[0];
	b = store->Sij[1];
	c = store->Sij[2];

	real disc = sqrt((a - c) * (a - c) + 4.0 * b * b);
	real r1 = ((a + c) + copysign(disc, (a - c))) / 2.0;
	real r2 = ((a + c) - copysign(disc, (a - c))) / 2.0;
	real r;
	real Fx, Fy;
	real phi;

	if ((a - c) == 0.0)
		phi = 0;
	else
		phi = atan(2 * b / (a - c));

	if (fabs(store->ux) > fabs(store->uy)) {
		r = r1;
		Fx = cos(phi) * copysign(1, -store->ux);
		Fy = sin(phi) * copysign(1, -store->ux);
	}
	else {
		r = r2;
		Fx = -sin(phi) * copysign(1, -store->uy);
		Fy = cos(phi) * copysign(1, -store->uy);
	}

	Fx *= ps->margForceRatio * r;
	Fy *= ps->margForceRatio * r;

	/*
	Fx = b/sqrt(b*b - (r-a)*(r-a));
	Fy = (r-a)/sqrt(b*b-(r-a)*(r-a));

	real Fmag = sqrt(Fx*Fx + Fy*Fy);

	if(fabs(Fmag) < FLT_EPSILON ) {
		Fx = 0.0;
		Fy = 0.0;
	}
	else {
		Fx *= rho/Fmag*ps->margForceRatio*abs(r);
		Fy *= rho/Fmag*ps->margForceRatio*abs(r);
	}

	*/

	//Fx *= 1e-22;
	//Fy *= 1e-22;

	// TODO: calculate proper forces!!!
/*
	Fx = cos(smax_ang)*smax_mag*100.0f;
	Fy = sin(smax_ang)*smax_mag*100.0f;

	real Fu = Fx*store->ux+Fy*store->uy;

    int iPop;
    for(iPop=0; iPop<9; iPop++) {
        fPop[iPop] *= (1-omega);
        fPop[iPop] += omega * computePSEquilibrium (iPop, rho, store->ux, store->uy) +
        		computeForce(iPop, rho, omega, store->ux, store->uy, Fx, Fy, Fu);
    }
*/



	// TODO: -sin: the minus is a stupid correction. check why is it necessary
	//Fx = rho*cos(smax_ang)*smax_mag*ps->margForceRatio;
	//Fy = rho*sin(smax_ang)*smax_mag*ps->margForceRatio;

	real Fu = Fx * store->ux + Fy * store->uy;

	real forcePop[9];
	int iPop = 0;
	for (iPop = 0; iPop < 9; iPop++)
		forcePop[iPop] = computeForce(iPop, rho, omega, store->ux, store->uy, Fx, Fy, Fu);


	real n_ux, n_uy;
	n_ux = store->ux + Fx / 2.0 / rho;
	n_uy = store->uy + Fy / 2.0 / rho;

	// Store couple field force instead of velocity field force
	store->fx = Fx;
	store->fy = Fy;

	real omega_min = (1 - omega);

	for (iPop = 0; iPop < 9; iPop++) {
		fPop[iPop] *= omega_min;
		fPop[iPop] += omega * computePSEquilibrium(iPop, rho, n_ux, n_uy) + forcePop[iPop];
	}
}


/* BGK Dynamics                                                    */
/*******************************************************************/

// compute local equilibrium from rho and u
real inline computeEquilibrium(int iPop, real rho, real ux, real uy, real uSqr) {
	real c_u = c[iPop][0] * ux + c[iPop][1] * uy;
	return rho * t[iPop] * (1. + 3. * c_u + 4.5 * c_u * c_u - 1.5 * uSqr);
}

// bgk collision term
void inline bgk(real *fPop, void *selfData, NodeStore *store) {

	fluidData *fp = (fluidData *) selfData;

	real rho, ux, uy, uSqr;
	real omega = fp->omega;

    real f_eq[9];

	computeMacros(fPop, &rho, &ux, &uy);

	// Store macroscopic velocity (for use with coupled fields)
	store->ux = ux;
	store->uy = uy;

	if (fp->computeStressTensor) {
		store->Sij[0] = 0;
		store->Sij[1] = 0;
		store->Sij[2] = 0;
	}

	uSqr = ux * ux + uy * uy;

    // Compute equilibrium and non-equilibrium part of the stress tensor
	int iPop;
	for (iPop = 0; iPop < 9; iPop++) {
		f_eq[iPop] = computeEquilibrium(iPop, rho, ux, uy, uSqr);

		if (fp->computeStressTensor) {
			real f_neq = fPop[iPop] - f_eq[iPop];
			store->Sij[0] += f_neq * c[iPop][0] * c[iPop][0];
			store->Sij[1] += f_neq * c[iPop][0] * c[iPop][1];
			store->Sij[2] += f_neq * c[iPop][1] * c[iPop][1];
		}
	}

    // Compute the strain rate tensor from the stress tensor
	if (fp->computeStressTensor) {
		//real sigmaFact =  - omega / (1./(3.0*omega) + 0.1666666);    //TODO: i'm not sure about this...
		real sigmaFact = -omega / (2. / 3. * rho);
		store->Sij[0] *= sigmaFact;
		store->Sij[1] *= sigmaFact;
		store->Sij[2] *= sigmaFact;
		store->mu = rho * (1.0 / (3.0 * omega) - 1.0 / 6.0);
	}

	// Recompute local omega
	if (fp->useNonNewt) {
		NonNewtParam *np = fp->np;

		real shearRate = sqrt(store->Sij[0] * store->Sij[0] +
				2 * store->Sij[1] * store->Sij[1] +
				store->Sij[2] * store->Sij[2]);

		real nu = np->nu_inf + (np->nu_0 - np->nu_inf) * pow(1 + pow(np->lambda * shearRate, np->a), (np->n - 1) / np->a);
		omega = 1.0 / (3.0 * nu + 0.5);
        fp->omega = omega;
	}

    // Apply Smagorinsky turbulence model
    if(fp->useSmagorinsky) {
        real tau0 = 1./omega;
        real sqrtStress = sqrt(store->Sij[0]*store->Sij[0]+2.0*store->Sij[1]*store->Sij[1]+store->Sij[2]*store->Sij[2]);
        real tau_t = 0.5*(sqrt(tau0*tau0 + 18.0*fp->smagorinskyC*fp->smagorinskyC*sqrtStress)-tau0);

        omega = 1.0/(tau0+tau_t);
    }

    real omega_min = 1.0 - omega;

    // Collide
    for(iPop=0; iPop<9; iPop++)
    {
        fPop[iPop] *= omega_min;
        fPop[iPop] += omega * f_eq[iPop];
    }
}

/*  MRT Dynamics                                                    */
/********************************************************************/

// compute local equilibrium from rho and u
void inline computeMRTEquilibrium(real *m_eq, real rho, real ux, real uy) {
	real uSqr = ux * ux + uy * uy;

	m_eq[0] = rho;                        // density
	m_eq[1] = rho * (-2 + 3*uSqr);        // energy
	m_eq[2] = rho * (1 - 3*uSqr);         // energy square
	m_eq[3] = rho * ux;                   // momentum in x-dir
	m_eq[4] = rho * (-ux);                // energy flux in x-dir
	m_eq[5] = rho * uy;                   // momentum in y-dir
	m_eq[6] = rho * (-uy);                // energy flux in y-dir
	m_eq[7] = rho * (ux * ux - uy * uy);  // diagonal component of stress tensor
	m_eq[8] = rho * ux * uy;              // off-diagonal component of stress tensor
}

void inline bgk2mrt(real *fPop, real *mPop) {
	int k, l;

	for (k = 0; k < 9; k++) {
		mPop[k] = 0.0;
		for (l = 0; l < 9; l++)
			mPop[k] += M[k][l] * fPop[l];
	}
}

void inline mrt2bgk(real *mPop, real *fPop) {
	int k, l;

	for (k = 0; k < 9; k++) {
		fPop[k] = 0.0;
		for (l = 0; l < 9; l++)
			fPop[k] += Minv[k][l] * mPop[l];
	}
}


void inline mrtCollision(real *m, real *mEq, void *selfData) {
	fluidData *fp = (fluidData *) selfData;

    //real s_q = (8*fp->omega-16.0)/(fp->omega-8.0);

	mEq[0] = 0;//mEq[0]=(m[0]-mEq[0])*S[0];
	mEq[1] = (m[1] - mEq[1]) * S[1];
	mEq[2] = (m[2] - mEq[2]) * S[2];
	mEq[3] = 0;//mEq[3]=(m[3]-mEq[3])*S[3];
	mEq[4] = (m[4] - mEq[4]) * S[4];//s_q;//S[4];
	mEq[5] = 0;//mEq[5]=(m[5]-mEq[5])*S[5];
	mEq[6] = (m[6] - mEq[6]) * S[6];//s_q;//S[6];
	mEq[7] = (m[7] - mEq[7]) * fp->omega;
	mEq[8] = (m[8] - mEq[8]) * fp->omega;
}

// mrt collision term - TODO: shear stress computation!!!
void inline mrt(real *fPop, void *selfData, NodeStore *store) {
	real rho, ux, uy;

	fluidData *fp = (fluidData *) selfData;

	real m_coll[9];
	real m[9];

    //real omega = fp->omega;

	computeMacros(fPop, &rho, &ux, &uy);

	store->ux = ux;
	store->uy = uy;

	if (fp->computeStressTensor) {
		store->Sij[0] = 0;
		store->Sij[1] = 0;
		store->Sij[2] = 0;
	}

	bgk2mrt(fPop, m); // m = M*fPop

	computeMRTEquilibrium(m_coll, rho, ux, uy); // m_coll = m_eq

	if (fp->computeStressTensor) {
		// change back to f_eq + use fPop
		real fEq[9];
		mrt2bgk(m_coll, fEq);

		int iPop;
		for (iPop = 0; iPop < 9; iPop++) {
			real f_neq = fPop[iPop] - fEq[iPop];
			store->Sij[0] += f_neq * c[iPop][0] * c[iPop][0];
			store->Sij[1] += f_neq * c[iPop][0] * c[iPop][1];
			store->Sij[2] += f_neq * c[iPop][1] * c[iPop][1];
		}

		//real sigmaFact =  - omega / (1./(3.0*omega) + 0.1666666);    //TODO: i'm not so sure about this...
		real sigmaFact = -fp->omega / (2. / 3. * rho);
		store->Sij[0] *= sigmaFact;
		store->Sij[1] *= sigmaFact;
		store->Sij[2] *= sigmaFact;
		store->mu = rho * (1.0 / (3.0 * fp->omega) - 1.0 / 6.0);
	}

    // Recompute local omega based on rheology
    if (fp->useNonNewt) {
        NonNewtParam *np = fp->np;

        real shearRate = sqrt(store->Sij[0] * store->Sij[0] +
                              2 * store->Sij[1] * store->Sij[1] +
                              store->Sij[2] * store->Sij[2]);

        real nu = np->nu_inf + (np->nu_0 - np->nu_inf) * pow(1 + pow(np->lambda * shearRate, np->a), (np->n - 1) / np->a);
        fp->omega = 1.0 / (3.0 * nu + 0.5);
    }

    // Apply Smagorinsky turbulence model
    if(fp->useSmagorinsky) {
        real tau0 = 1./fp->omega;
        real sqrtStress = sqrt(store->Sij[0]*store->Sij[0]+2.0*store->Sij[1]*store->Sij[1]+store->Sij[2]*store->Sij[2]);
        real tau_t = 0.5*(sqrt(tau0*tau0 + 18.0*fp->smagorinskyC*fp->smagorinskyC*sqrtStress)-tau0);

        fp->omega = 1.0/(tau0+tau_t);
    }

    mrtCollision(m, m_coll, selfData);  // m_coll = S*(m-m_eq)

	int i;
	for (i = 0; i < 9; i++)                    //m = m - m_coll
		m[i] -= m_coll[i];

	mrt2bgk(m, fPop); //fPop = M^-1 * m


}

void inline mrtForced(real *fPop, void *selfData, NodeStore *store) {

    porousData *pd = (porousData *) selfData;

    real omega = pd->fp.omega;

    real rho, ux, uy, uSqr, Fx, Fy;

    real forcePop[9];
    real f_eq[9];

    computeMacros(fPop, &rho, &ux, &uy);

    // Calculate forcing term
    porousForcing(&Fx, &Fy, pd->c0, pd->c1, ux, uy);

    real Fu = Fx * ux + Fy * uy;

    int iPop;
    for (iPop = 0; iPop < 9; iPop++)
        forcePop[iPop] = computeForce(iPop, rho, omega, ux, uy, Fx, Fy, Fu);

    // Modify macroscopic velocity
    ux += Fx / 2.0 / rho;
    uy += Fy / 2.0 / rho;

    // Store macroscopic values
    store->ux = ux;
    store->uy = uy;

    store->fx = Fx;
    store->fy = Fy;


	fluidData *fp = &pd->fp;

	real m_coll[9];
	real m[9];

	//real omega = fp->omega;

	computeMacros(fPop, &rho, &ux, &uy);

	store->ux = ux;
	store->uy = uy;

	if (fp->computeStressTensor) {
		store->Sij[0] = 0;
		store->Sij[1] = 0;
		store->Sij[2] = 0;
	}

	bgk2mrt(fPop, m); // m = M*fPop

	computeMRTEquilibrium(m_coll, rho, ux, uy); // m_coll = m_eq

	if (fp->computeStressTensor) {
		// change back to f_eq + use fPop
		real fEq[9];
		mrt2bgk(m_coll, fEq);

		int iPop;
		for (iPop = 0; iPop < 9; iPop++) {
			real f_neq = fPop[iPop] - fEq[iPop];
			store->Sij[0] += f_neq * c[iPop][0] * c[iPop][0];
			store->Sij[1] += f_neq * c[iPop][0] * c[iPop][1];
			store->Sij[2] += f_neq * c[iPop][1] * c[iPop][1];
		}

		//real sigmaFact =  - omega / (1./(3.0*omega) + 0.1666666);    //TODO: i'm not so sure about this...
		real sigmaFact = -fp->omega / (2. / 3. * rho);
		store->Sij[0] *= sigmaFact;
		store->Sij[1] *= sigmaFact;
		store->Sij[2] *= sigmaFact;
		store->mu = rho * (1.0 / (3.0 * fp->omega) - 1.0 / 6.0);
	}

	// Recompute local omega based on rheology
	if (fp->useNonNewt) {
		NonNewtParam *np = fp->np;

		real shearRate = sqrt(store->Sij[0] * store->Sij[0] +
							  2 * store->Sij[1] * store->Sij[1] +
							  store->Sij[2] * store->Sij[2]);

		real nu = np->nu_inf + (np->nu_0 - np->nu_inf) * pow(1 + pow(np->lambda * shearRate, np->a), (np->n - 1) / np->a);
		fp->omega = 1.0 / (3.0 * nu + 0.5);
	}

	// Apply Smagorinsky turbulence model
	if(fp->useSmagorinsky) {
		real tau0 = 1./fp->omega;
		real sqrtStress = sqrt(store->Sij[0]*store->Sij[0]+2.0*store->Sij[1]*store->Sij[1]+store->Sij[2]*store->Sij[2]);
		real tau_t = 0.5*(sqrt(tau0*tau0 + 18.0*fp->smagorinskyC*fp->smagorinskyC*sqrtStress)-tau0);

		fp->omega = 1.0/(tau0+tau_t);
	}

	mrtCollision(m, m_coll, selfData);  // m_coll = S*(m-m_eq)

	int i;
	for (i = 0; i < 9; i++)                    //m = m - m_coll
		m[i] -= m_coll[i];

	mrt2bgk(m, fPop); //fPop = M^-1 * m

    // Add forcing distribution components
    for(iPop=0; iPop < 9; iPop++)
        fPop[iPop] += forcePop[iPop];

}


void inline ageingDyn(real *fPop, void *selfData, NodeStore *store)
{

    ageData *ad = (ageData*)selfData;
    real *fluidPop = ad->fluidPop;

    real fluidRho;
    computeRho(fluidPop, &fluidRho);

    real age;
    computeRho(fPop, &age);
    age /= fluidRho;

	int i;
    for(i = 0; i < 9; i++)
        fPop[i] = (age + ad->dt) * fluidPop[i];

}
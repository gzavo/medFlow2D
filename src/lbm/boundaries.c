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


#include <stdlib.h>
#include <stdio.h>
#include "boundaries.h"
#include "lb.h"

/* D2Q9 lattice constants                                        */
/*****************************************************************/
/*
  // opposite directions, for bounce back implementation
static const int oppositeOf[9] = { 0, 5, 6, 7, 8, 1, 2, 3, 4 };

  // lattice weights
static const real t[9] = { 4./9., 1./9., 1./36., 1./9., 1./36.,
                             1./9., 1./36., 1./9., 1./36. };
  // lattice velocities
static const int c[9][2] = {
    { 0, 0},
    { 1, 0}, { 1, 1}, {0,  1}, {-1, 1},
    {-1, 0}, {-1,-1}, {0, -1}, {1, -1}
};
*/


/* No-dynamics                                                   */
/*****************************************************************/

void inline nodynamics(real* fPop, void* selfData, NodeStore* store) {

}


/* Bounce back                                                   */
/*****************************************************************/
// TODO: optimise - loop unroll
void inline bounceBack(real* fPop, void* selfData, NodeStore* store) {
    real fTmp[9];
    int iPop;
    for (iPop=1; iPop<9; iPop++) {
        fTmp[iPop] = fPop[oppositeOf[iPop]];
    }
    for (iPop=1; iPop<9; iPop++) {
        fPop[iPop] = fTmp[iPop];
    }
}

/* Helper functions: compute rho from u and vice versa           */
/*****************************************************************/

  /* Compute density on wall from bulk information on
     upper boundary. */
inline static real topRho(real *fPop, real uy) {
    return (fPop[0] + fPop[1] + fPop[5] + 2*(fPop[2]+fPop[3]+fPop[4])) / (1.-uy);
}

  /* Compute uy on wall from bulk information on
     upper boundary. */
inline static real topU(real *fPop, real rho) {
    return 1. - (fPop[0] + fPop[5] + fPop[1] + 2*(fPop[4]+fPop[3]+fPop[2])) / rho;
}

  /* Compute density on wall from bulk information on
     lower boundary. */
inline static real bottomRho(real *fPop, real uy) {
    return (fPop[0] + fPop[1] + fPop[5] + 2*(fPop[6]+fPop[7]+fPop[8])) / (1+uy);
}

  /* Compute uy on wall from bulk information on
     lower boundary. */
inline static real bottomU(real *fPop, real rho) {
    return -1. + (fPop[0] + fPop[1] + fPop[5] + 2*(fPop[6]+fPop[7]+fPop[8])) / rho;
}

  /* Compute density on wall from bulk information on
     right boundary. */
inline static real rightRho(real* fPop, real ux) {
    return (fPop[0] + fPop[3] + fPop[7] + 2*(fPop[1]+fPop[2]+fPop[8])) / (1+ux);
}

  /* Compute ux on wall from bulk information on
     right boundary. */
inline static real rightU(real* fPop, real rho) {
    return -1. + (fPop[0] + fPop[3] + fPop[7] + 2*(fPop[1]+fPop[2]+fPop[8])) / rho;
}

  /* Compute density on wall from bulk information on
     left boundary. */
inline static real leftRho(real* fPop, real ux) {
    return (fPop[0] + fPop[3] + fPop[7] + 2*(fPop[4]+fPop[5]+fPop[6])) / (1.-ux);
}

  /* Compute ux on wall from bulk information on
     left boundary. */
inline static real leftU(real* fPop, real rho) {
    return 1. - (fPop[0] + fPop[3] + fPop[7] + 2*(fPop[4]+fPop[5]+fPop[6])) / rho;
} 

/* Zou/He helper functions and boundary implementations          */
/*****************************************************************/

inline static void completeTop(real *fPop, real ux, real uy, real rho)
{
    fPop[6] = rho*uy/6 - rho*ux/2. + fPop[2] + 0.5 * (fPop[1]-fPop[5]);
    fPop[7] = fPop[3] + 2./3.*rho*uy;
    fPop[8] = rho*uy/6 + rho*ux/2. + fPop[4] + 0.5 *(fPop[5]-fPop[1]);

}

inline static void completeBottom(real *fPop, real ux, real uy, real rho)
{
    fPop[2] = - rho*uy/6 + rho*ux/2. + fPop[6] + 0.5 *(fPop[5]-fPop[1]);
    fPop[3] = fPop[7] - 2./3.*rho*uy;
    fPop[4] = - rho*uy/6 - rho*ux/2. + fPop[8] + 0.5 *(fPop[1]-fPop[5]);
}

inline static void completeRight(real* fPop, real ux, real uy, real rho)
{
    fPop[4] = - rho*ux/6 - rho*uy/2. + fPop[8] + 0.5 *(fPop[7]-fPop[3]);
    fPop[5] = fPop[1] - 2./3.*rho*ux;
    fPop[6] = - rho*ux/6 - rho*uy/2. + fPop[2] + 0.5 *(fPop[3]-fPop[7]);
}

inline static void completeLeft(real* fPop, real ux, real uy, real rho)
{
    fPop[1] = fPop[5] + 2./3.*rho*ux;
    fPop[2] = rho*ux/6. - rho*uy/2. + fPop[6] + 0.5 *(fPop[7]-fPop[3]);
    fPop[8] = rho*ux/6. + rho*uy/2. + fPop[4] + 0.5 *(fPop[3]-fPop[7]);
}
 

void inline topZouHe(real *fPop, void *selfData, NodeStore *store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real rho = topRho(fPop, data->uy);
    completeTop(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline bottomZouHe(real *fPop, void *selfData, NodeStore *store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real rho = bottomRho(fPop, data->uy);
    completeBottom(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline leftZouHe(real* fPop, void* selfData, NodeStore* store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real rho = leftRho(fPop, data->ux);
    completeLeft(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline rightZouHe(real* fPop, void* selfData, NodeStore* store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real rho = rightRho(fPop, data->ux);
    completeRight(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline topPressureZouHe(real *fPop, void *selfData, NodeStore *store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real uy = topU(fPop, data->rho);
    completeTop(fPop, data->uPar, uy, data->rho);
    data->bulkDynamics->dynamicsFun ( fPop, data->bulkDynamics->selfData, store );
}

void inline bottomPressureZouHe(real *fPop, void *selfData, NodeStore *store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real uy = bottomU(fPop, data->rho);
    completeBottom(fPop, data->uPar, uy, data->rho);
    data->bulkDynamics->dynamicsFun ( fPop, data->bulkDynamics->selfData, store );
}

void inline leftPressureZouHe(real* fPop, void* selfData, NodeStore* store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real ux = leftU(fPop, data->rho);
    completeLeft(fPop, ux, data->uPar, data->rho);
    data->bulkDynamics->dynamicsFun ( fPop, data->bulkDynamics->selfData, store );
}

void inline rightPressureZouHe(real* fPop, void* selfData, NodeStore* store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real ux = rightU(fPop, data->rho);
    completeRight(fPop, ux, data->uPar, data->rho);
    data->bulkDynamics->dynamicsFun ( fPop, data->bulkDynamics->selfData, store );
} 


/* Regularized helper functions and boundary implmenetations     */
/*****************************************************************/

inline static void splitEqNeq (real* f, real* fEq, real* fNeq, real rho, real ux, real uy)
{
    real uSqr = ux*ux+uy*uy;
    int iPop;
    for (iPop=0; iPop<9; iPop++) {
        fEq[iPop]  = computeEquilibrium(iPop, rho, ux, uy, uSqr);
        fNeq[iPop] = f[iPop]-fEq[iPop];
    }
}

inline static void regularizedF ( real* f, real* fEq, real *fNeq )
{
    real neqPixx=0, neqPiyy=0, neqPixy=0;

    int iPop;
    for (iPop = 0; iPop < 9; iPop++)
    {
        neqPixx += fNeq[iPop] * c[iPop][0] * c[iPop][0];
        neqPixy += fNeq[iPop] * c[iPop][0] * c[iPop][1];
        neqPiyy += fNeq[iPop] * c[iPop][1] * c[iPop][1];
    }

    for (iPop=0; iPop<9; iPop++) {
        f[iPop] = fEq[iPop] + 9./2. * t[iPop] *
            ( (c[iPop][0]*c[iPop][0]-1./3.)*neqPixx +
              (c[iPop][1]*c[iPop][1]-1./3.)*neqPiyy +
              2.*c[iPop][0]*c[iPop][1]*neqPixy );
    }
}

void inline topRegularized(real *fPop, void *selfData, NodeStore *store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real fEq[9], fNeq[9];

    real rho = topRho(fPop, data->uy);

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);

    fNeq[6] = fNeq[2];
    fNeq[7] = fNeq[3];
    fNeq[8] = fNeq[4];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline bottomRegularized(real *fPop, void *selfData, NodeStore *store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real fEq[9], fNeq[9];

    real rho = bottomRho(fPop, data->uy);

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);

    fNeq[2] = fNeq[6];
    fNeq[3] = fNeq[7];
    fNeq[4] = fNeq[8];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline leftRegularized(real* fPop, void* selfData, NodeStore* store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real fEq[9], fNeq[9];

    real rho = leftRho(fPop, data->ux);

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);

    fNeq[2] = fNeq[4];
    fNeq[1] = fNeq[5];
    fNeq[8] = fNeq[6];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline rightRegularized(real* fPop, void* selfData, NodeStore* store) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    real fEq[9], fNeq[9];

    real rho = rightRho(fPop, data->ux);

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);

    fNeq[4] = fNeq[2];
    fNeq[5] = fNeq[1];
    fNeq[6] = fNeq[8];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline topPressureRegularized(real *fPop, void *selfData, NodeStore *store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real fEq[9], fNeq[9];

    real uy  = topU(fPop, data->rho);

    splitEqNeq(fPop, fEq, fNeq, data->rho, data->uPar, uy);

    fNeq[6] = fNeq[2];
    fNeq[7] = fNeq[3];
    fNeq[8] = fNeq[4];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline bottomPressureRegularized(real *fPop, void *selfData, NodeStore *store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real fEq[9], fNeq[9];

    real uy  = bottomU(fPop, data->rho);

    splitEqNeq(fPop, fEq, fNeq, data->rho, data->uPar, uy);

    fNeq[2] = fNeq[6];
    fNeq[3] = fNeq[7];
    fNeq[4] = fNeq[8];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

void inline leftPressureRegularized(real* fPop, void* selfData, NodeStore* store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real fEq[9], fNeq[9];

    real ux = leftU(fPop, data->rho);

    splitEqNeq(fPop, fEq, fNeq, data->rho, ux, data->uPar);

    fNeq[2] = fNeq[4];
    fNeq[1] = fNeq[5];
    fNeq[8] = fNeq[6];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store);
}

void inline rightPressureRegularized(real* fPop, void* selfData, NodeStore* store) {
    PressureBCData* data = (PressureBCData*) selfData;
    real fEq[9], fNeq[9];

    real ux = rightU(fPop, data->rho);

    splitEqNeq(fPop, fEq, fNeq, data->rho, ux, data->uPar);

    fNeq[4] = fNeq[2];
    fNeq[5] = fNeq[1];
    fNeq[6] = fNeq[8];

    regularizedF(fPop, fEq, fNeq);

    data->bulkDynamics->dynamicsFun (fPop, data->bulkDynamics->selfData, store );
}

/******************** Passive scalar boundaries ***********/

void inline topPressurePS(real *fPop, void *selfData, NodeStore *store)
{
	PressureBCData* data = (PressureBCData*) selfData;

	real rhoResid = (data->rho - (fPop[0]+fPop[1]+fPop[2]+fPop[3]+fPop[4]+fPop[5]))/(t[6]+t[7]+t[8]);

	fPop[6] = t[6]*rhoResid;
	fPop[7] = t[7]*rhoResid;
	fPop[8] = t[8]*rhoResid;

	data->bulkDynamics->dynamicsFun(fPop, data->bulkDynamics->selfData, store);
}

void inline bottomPressurePS(real *fPop, void *selfData, NodeStore *store)
{
	PressureBCData* data = (PressureBCData*) selfData;

	real rhoResid = (data->rho - (fPop[0]+fPop[1]+fPop[5]+fPop[6]+fPop[7]+fPop[8]))/(t[2]+t[3]+t[4]);

	fPop[2] = t[2]*rhoResid;
	fPop[3] = t[3]*rhoResid;
	fPop[4] = t[4]*rhoResid;

	data->bulkDynamics->dynamicsFun(fPop, data->bulkDynamics->selfData, store);

}

void inline leftPressurePS(real* fPop, void* selfData, NodeStore* store)
{
	PressureBCData* data = (PressureBCData*) selfData;

	real rhoResid = (data->rho - (fPop[0]+fPop[3]+fPop[4]+fPop[5]+fPop[6]+fPop[7]))/(t[1]+t[2]+t[8]);

	fPop[1] = t[1]*rhoResid;
	fPop[2] = t[2]*rhoResid;
	fPop[8] = t[8]*rhoResid;

	data->bulkDynamics->dynamicsFun(fPop, data->bulkDynamics->selfData, store);

}

void inline rightPressurePS(real* fPop, void* selfData, NodeStore* store)
{
	PressureBCData* data = (PressureBCData*) selfData;

	real rhoResid = (data->rho - (fPop[0]+fPop[1]+fPop[2]+fPop[3]+fPop[7]+fPop[8]))/(t[4]+t[5]+t[6]);

	fPop[4] = t[4]*rhoResid;
	fPop[5] = t[5]*rhoResid;
	fPop[6] = t[6]*rhoResid;

	data->bulkDynamics->dynamicsFun(fPop, data->bulkDynamics->selfData, store);

}


void inline zeroAge(real* fPop, void* selfData, NodeStore* store)
{

    int i;
    for(i=0; i < 9; i++)
        fPop[i] = 0;

}

void inline leftAgeOutlet(real* fPop, void* selfData, NodeStore* store)
{
    ageData *ad = (ageData*)selfData;
    real *fluidPop = ad->fluidPop;

    real rhoPartial = fluidPop[0] + fluidPop[3] + fluidPop[4] + fluidPop[5] + fluidPop[6] + fluidPop[7];
    real agePartial = fPop[0] + fPop[3] + fPop[4] + fPop[5] + fPop[6] + fPop[7];

    real age = agePartial / rhoPartial;

    int i;
    for(i = 0; i < 9; i++)
        fPop[i] = (age + ad->dt) * fluidPop[i];
}

void inline rightAgeOutlet(real* fPop, void* selfData, NodeStore* store)
{
    ageData *ad = (ageData*)selfData;
    real *fluidPop = ad->fluidPop;

    real rhoPartial = fluidPop[0] + fluidPop[1] + fluidPop[2] + fluidPop[3] + fluidPop[7] + fluidPop[8];
    real agePartial = fPop[0] + fPop[1] + fPop[2] + fPop[3] + fPop[7] + fPop[8];

    real age = agePartial / rhoPartial;

    int i;
    for(i = 0; i < 9; i++)
        fPop[i] = (age + ad->dt) * fluidPop[i];
}

void inline topAgeOutlet(real* fPop, void* selfData, NodeStore* store)
{
    ageData *ad = (ageData*)selfData;
    real *fluidPop = ad->fluidPop;

    real rhoPartial = fluidPop[0] + fluidPop[1] + fluidPop[5] + fluidPop[6] + fluidPop[7] + fluidPop[8];
    real agePartial = fPop[0] + fPop[1] + fPop[5] + fPop[6] + fPop[7] + fPop[8];

    real age = agePartial / rhoPartial;

    int i;
    for(i = 0; i < 9; i++)
        fPop[i] = (age + ad->dt) * fluidPop[i];
}

void inline bottomAgeOutlet(real* fPop, void* selfData, NodeStore* store)
{
    ageData *ad = (ageData*)selfData;
    real *fluidPop = ad->fluidPop;

    real rhoPartial = fluidPop[0] + fluidPop[1] + fluidPop[2] + fluidPop[3] + fluidPop[4] + fluidPop[5];
    real agePartial = fPop[0] + fPop[1] + fPop[2] + fPop[3] + fPop[4] + fPop[5];

    real age = agePartial / rhoPartial;

    int i;
    for(i = 0; i < 9; i++)
        fPop[i] = (age + ad->dt) * fluidPop[i];
}

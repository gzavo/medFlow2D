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



#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "lb.h"

typedef enum { UNSET = -1, LEFT, TOP, RIGHT, BOTTOM } domainSide;


/* struct VelocityBCData                                         */
/*****************************************************************/
/*   This struct is used as a parameter to the velocity boundary
 *   conditions. It takes place of the variable "selfData" in the
 *   struct "Dynamics". It contains the value of the Dirichlet
 *   boundary, and a pointer to the underlying bulk dynamics, e.g.
 *   bgk.
 */
typedef struct {
    Dynamics* bulkDynamics;
    real    ux, uy;
} VelocityBCData;

/* struct PressureBCData                                         */
/*****************************************************************/
/*   This struct is used as a parameter to the pressure boundary
 *   conditions. It takes place of the variable "selfData" in the
 *   struct "Dynamics". It contains the value of the velocity
 *   tangential to the boundary, the value of the prescribed
 *   pressure, and a pointer to the underlying bulk dynamics, e.g.
 *   bgk.
 */
typedef struct {
    Dynamics* bulkDynamics;
    real    rho, uPar;
} PressureBCData;

/* All the implemented boundaries...                             */
/*****************************************************************/

  // no-dynamics for unused cells
void nodynamics(real* fPop, void* selfData, NodeStore* store);

  // bounce back no-slip condition. "selfData" is a null pointer,
  //   given that the function need no parameter.
void bounceBack(real* fPop, void* selfData, NodeStore* store);

  // ZouHe velocity boundaries on upper, lower, left and right
  //   boundaries
void topZouHe(real *fPop, void *selfData, NodeStore *store);
void bottomZouHe(real *fPop, void *selfData, NodeStore *store);
void leftZouHe (real* fPop, void* selfData, NodeStore* store);
void rightZouHe(real* fPop, void* selfData, NodeStore* store);

  // ZouHe pressure boundaries on upper, lower, left and right
  //   boundaries
void topPressureZouHe(real *fPop, void *selfData, NodeStore *store);
void bottomPressureZouHe(real *fPop, void *selfData, NodeStore *store);
void leftPressureZouHe (real* fPop, void* selfData, NodeStore* store);
void rightPressureZouHe(real* fPop, void* selfData, NodeStore* store);
 
  // Regularized velocity boundaries on upper, lower, left and
  // right boundaries
void topRegularized(real *fPop, void *selfData, NodeStore *store);
void bottomRegularized(real *fPop, void *selfData, NodeStore *store);
void leftRegularized (real* fPop, void* selfData, NodeStore* store);
void rightRegularized(real* fPop, void* selfData, NodeStore* store);

  // Regularized pressure boundaries on upper, lower, left and
  // right boundaries
void topPressureRegularized(real *fPop, void *selfData, NodeStore *store);
void bottomPressureRegularized(real *fPop, void *selfData, NodeStore *store);
void leftPressureRegularized (real* fPop, void* selfData, NodeStore* store);
void rightPressureRegularized(real* fPop, void* selfData, NodeStore* store);

  // Passive scalar constant density conditions
void topPressurePS(real *fPop, void *selfData, NodeStore *store);
void bottomPressurePS(real *fPop, void *selfData, NodeStore *store);
void leftPressurePS(real* fPop, void* selfData, NodeStore* store);
void rightPressurePS(real* fPop, void* selfData, NodeStore* store);


  // Water age boundaries
void zeroAge(real* fPop, void* selfData, NodeStore* store);

void leftAgeOutlet(real* fPop, void* selfData, NodeStore* store);

void rightAgeOutlet(real* fPop, void* selfData, NodeStore* store);

void topAgeOutlet(real* fPop, void* selfData, NodeStore* store);

void bottomAgeOutlet(real* fPop, void* selfData, NodeStore* store);

#endif

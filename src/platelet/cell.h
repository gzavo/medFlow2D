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


#ifndef CELL_H
#define CELL_H

#include "../common.h"
#include "../utils/arraylist.h"
#include "../utils/queue.h"
#include "../lbm/lb.h"

typedef struct {
	real shearStress;
	real rho;
} CellHistoryData;

typedef struct {
	int x;
	int y;
	Queue *queue;
	int queueLength;
	real shearStressSum;
	real rhoSum;
	real rhoADP;
} CellHistory;

int cellIndexEqual(ArrayListValue value1, ArrayListValue value2);
int cellIndexCompare(ArrayListValue value1, ArrayListValue value2);

inline real coagProbabilityFunction(real rhoADP, real rhoPlateletSum, real shearStressSum);
void setRhoADP(ArrayList *a, int iX, int iY, real newRhoADP);

int getMaxHistory();
void setMaxHistory(int maxHist);
void setRhoADPDecrement(real prob);
void setDefaultRhoADP(real defProb);
real getDefaultRhoADP();
void setCoagCoeff(real newCoagCoeff);
real getCoagCoeff();

ArrayList *buildNearWallCells(Simulation *sim, Node** lattice);
void updateNearWallListHistory(Simulation *sim, Node **lattice, ArrayList *a);
void updateNearWallCells(Simulation *sim, ArrayList *a);

int modifyNearWallCellsPS(Simulation *sim, int width, void   (*dynamicsFun) (real* fPop, void* selfData, NodeStore* store));

CellHistory *createCellHistory();
void freeCellHistory(CellHistory* ch);
void freeNearWallCells(ArrayList* a);

void saveNearWallData(char* fName, ArrayList *a);

#endif



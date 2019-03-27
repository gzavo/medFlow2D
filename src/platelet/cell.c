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
#include <math.h>
#include <stdio.h>
#include "cell.h"
#include "../lbm/boundaries.h"
#include "../lbm/lb.h"


int maxHistory = 10;                // Maximum history to follow per near wall node
real defaultRhoADP = 3.0e-7;		// Default concentration of ADP in blood [mg/m^3]
real coagCoeff = 1.0;               // Balance of coagulation probability arguments

// TODO: This should be scaled with physical length
real rhoADPDecrement = 0.9;     // Activated platelets release ADP, but somewhat less -> decay in space

// Inlined definitions
extern inline real coagProbabilityFunction(real rhoADP, real rhoPlateletSum, real shearStressSum);

int getMaxHistory()
{
    return maxHistory;
}

void setMaxHistory(int maxHist)
{
	maxHistory = maxHist;
}

void setRhoADPDecrement(real prob)
{
	rhoADPDecrement = prob;
}

void setCoagCoeff(real newCoagCoeff)
{
	coagCoeff = newCoagCoeff;
}

real getCoagCoeff() {
	return coagCoeff;
}


void setDefaultRhoADP(real defProb)
{
	defaultRhoADP = defProb;
}

real getDefaultRhoADP() {
	return defaultRhoADP;
}

void setRhoADP(ArrayList *a, int iX, int iY, real newRhoADP) {
	
	CellHistory ch_tmp;
	ch_tmp.x = iX; ch_tmp.y = iY;
	
	int idx = arraylist_index_of(a, cellIndexEqual, (void*)&ch_tmp);
	if(idx<0)
		printf(" -- WARNING -- NWC cell not found!");
	else {
		CellHistory *ch = (CellHistory*)a->data[idx];
		ch->rhoADP = newRhoADP;
	}
}

void saveNearWallData(char* fName, ArrayList *a)
{
	FILE* oFile = fopen(fName, "w");
	
    fprintf(oFile, "# Number of NWCs: %d\n", a->length);

    fprintf(oFile, "# X Y coagProb rhoADP rhoPlateletAvg shearMaxAvg\n");

	int i;
	for(i = 0; i < a->length; ++i)
	{
		CellHistory *ch = (CellHistory*)a->data[i];
		real totalProbability = coagProbabilityFunction(ch->rhoADP, ch->rhoSum, ch->shearStressSum);
		int iX = ch->x;
		int iY = ch->y;
        if (ch->queueLength >= maxHistory) 
            fprintf(oFile, "%d %d %f %lf %f %f\n", iX, iY, totalProbability, ch->rhoADP, ch->rhoSum, ch->shearStressSum);
        else
            fprintf(oFile, "#### %d %d %lf %f %f %f\n", iX, iY, totalProbability, ch->rhoADP, ch->rhoSum, ch->shearStressSum);
		
	}
	
	fclose(oFile);	
}

inline bool hasNeighbor(Simulation *sim, int iX, int iY, NodeType neighbor)
{
    int dX, dY;
	for(dX = -1; dX <= 1; dX++ )
        for(dY = -1; dY <= 1; dY++)
        {
            if (sim->geometry[iX+dX][iY+dY]==neighbor) 
                return true;
        }

    return false;
}

inline real coagProbabilityFunction(real rhoADP, real rhoPlateletSum, real shearStressSum)
{
    return coagCoeff * (rhoADP * rhoPlateletSum / shearStressSum);
}

void updateNearWallCells(Simulation *sim, ArrayList *a)
{
	int i;
	for(i = 0; i < a->length; ++i)
	{
		CellHistory *ch = (CellHistory*)a->data[i];
		
		if(ch->queueLength < maxHistory)
			continue;

		// The average of the history
		real totalProbability = coagProbabilityFunction(ch->rhoADP, ch->rhoSum, ch->shearStressSum);
		
		if(totalProbability > 1.0) {
			int iX = ch->x;
			int iY = ch->y;
			
			// set dynamics to bounce back
			Dynamics *dyn   = (Dynamics*)calloc(1, sizeof(Dynamics));
			Dynamics *dynPS = (Dynamics*)calloc(1, sizeof(Dynamics));

			dyn->dynamicsFun = &bounceBack;
			dynPS->dynamicsFun = &bounceBack;

			// Free up previous dynamics
			free(sim->latticeA[iX][iY].dynamics);
			free(sim->psLatticeA[iX][iY].dynamics);
			if(sim->doubleBuffered) {
				free(sim->latticeB[iX][iY].dynamics);
				free(sim->psLatticeB[iX][iY].dynamics);
			}
			// TODO: selfdata is leaking memory

			// Set new ones
			setDynamics(sim, iX,iY, dyn);
			setPSDynamics(sim, iX, iY, dynPS);

			setGeometry(sim, iX,iY, COAGWALL);
			//sim->lattice[iX][iY].dynamics = dyn;
			//sim->tmpLattice[iX][iY].dynamics = dynPS;

			//sim->geometry[iX][iY] = COAGWALL;
			
            real parentRhoADP = ch->rhoADP;

			//remove from NWC list
			freeCellHistory(ch);
			arraylist_remove(a, i);
			
			// search for possible new NWCs
			int dX, dY;
			for(dX = -1; dX <= 1; dX++ )
				for(dY = -1; dY <= 1; dY++)
					if(sim->geometry[iX+dX][iY+dY] == FLUID)
					{
						CellHistory* ch_temp = createCellHistory();
						ch_temp->x = iX+dX; ch_temp->y = iY+dY;
						
                        //if (hasNeighbor(sim,iX+dX,iY+dY,WALL)) 
                        //    ch_temp->rhoADP = defaultRhoADP;
                        //else
                        //    ch_temp->rhoADP = defaultRhoADP*2;	// Platlets thether more easily to eachother

                        ch_temp->rhoADP = defaultRhoADP + parentRhoADP * rhoADPDecrement;
                        //ch_temp->rhoADP = parentRhoADP*rhoADPDecrement;
						
						arraylist_append(a, (void*)ch_temp);
						
						setGeometry(sim, iX+dX, iY+dY, NWFLUID);
					}	
		}
				
	}
	
}

void updateNearWallListHistory(Simulation *sim, Node **lattice, ArrayList *a)
{
	int i;
	for(i = 0; i < a->length; ++i)
	{
		CellHistory *ch = (CellHistory*)a->data[i];
		int iX = ch->x;
		int iY = ch->y;
		
		real s[3];
		real mu2 = 2.0 * sim->store[iX][iY].mu;
		s[0] = mu2 * sim->store[iX][iY].Sij[0];
		s[1] = mu2 * sim->store[iX][iY].Sij[1];
		s[2] = mu2 * sim->store[iX][iY].Sij[2];

		// Maximal shear stress and direction from Mohr's diagram
		real sqrDist = (s[0]-s[2])/2.0;
		real r = sqrt(sqrDist*sqrDist + s[1]*s[1]);
		
		real *f = lattice[iX][iY].fPop;
		real rho = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8];
		
		// Insert new data to history
		CellHistoryData *chd = (CellHistoryData*)calloc(1, sizeof(CellHistoryData));
		
		chd->shearStress = r;
		chd->rho = rho;
		
		queue_push_head(ch->queue, (void *)chd);
		ch->queueLength++;
		
		ch->shearStressSum += r;
		ch->rhoSum += rho;
		
		// if history is too lengthy, pop last element
		if(ch->queueLength > maxHistory)
		{
			CellHistoryData* ch_temp = (CellHistoryData*)queue_pop_tail(ch->queue);
			ch->shearStressSum -= ch_temp->shearStress;
			ch->rhoSum -= ch_temp->rho;
			ch->queueLength--;
			free(ch_temp);
		}
	}
	
}

  //Compare function
int cellIndexCompare(ArrayListValue value1, ArrayListValue value2)
{
  	CellHistory *v1 = (CellHistory*)value1;
  	CellHistory *v2 = (CellHistory*)value2;
  
  	if(v1->x < v2->x)
  		return -1;
  	else if (v1->x > v2->x)
  		return 1;
  	else {
  		if(v1->y < v2->y)
  			return -1;
  		else if (v1->y > v2->y)
  			return 1;
  		else
  			return 0;
  	}
}
  
  	// Equal function
int cellIndexEqual(ArrayListValue value1, ArrayListValue value2)
{
	CellHistory *v1 = (CellHistory*)value1;
  	CellHistory *v2 = (CellHistory*)value2;
  
  	//if((v1->x == v2->x) && (v1->y == v2->y))
	//	return 0;
  
	//return -1;
	return (v1->x == v2->x) && (v1->y == v2->y);
}

CellHistory *createCellHistory()
{
	CellHistory* ch = (CellHistory*)calloc(1, sizeof(CellHistory));
	ch->queue = queue_new();
	return ch;
}

void freeCellHistory(CellHistory* ch)
{
	while (!queue_is_empty(ch->queue)) {
		free(queue_pop_head(ch->queue));
	}
	queue_free(ch->queue);
	free(ch);
	ch = NULL;
}

void freeNearWallCells(ArrayList* a)
{
	int i;
	for(i=0; i < a->length; i++)
	{		
		freeCellHistory(a->data[i]);		
	}
	arraylist_free(a);
}

int modifyNearWallCellsPS(Simulation *sim, int width, void   (*dynamicsFun) (real* fPop, void* selfData, NodeStore* store))
{

	int iX, iY;

	// NOTE: Nothing can happen on the boundary!
	for(iX = 2; iX < sim->lx; iX++ )
		for(iY = 2; iY < sim->ly; iY++)
		{
			if(sim->geometry[iX][iY] == FLUID)
			{
				int dX, dY;
				bool isNWC = false;

				for(dX = -width; dX <= width; dX++ )
				{
					for(dY = -width; dY <= width; dY++)
					{
						if(sim->geometry[iX+dX][iY+dY] == WALL || sim->geometry[iX+dX][iY+dY] == COAGWALL)
						{
							isNWC = true;
							//dY = 2;
							break;
						}
					}
					if(isNWC)
						//dX = 2;
						break;
				}
				if(isNWC) {
					sim->psLatticeA[iX][iY].dynamics->dynamicsFun = dynamicsFun;
					if(sim->doubleBuffered)
						sim->psLatticeB[iX][iY].dynamics->dynamicsFun = dynamicsFun;
				}
			}
		}

	return 0;
}


ArrayList *buildNearWallCells(Simulation *sim, Node** lattice) 
{
	ArrayList *nwc;
	
	nwc = arraylist_new((sim->lx+sim->ly)*2);		// sensible bet

	int iX, iY;
	
	// NOTE: Nothing can happen on the boundary!
	for(iX = 2; iX < sim->lx; iX++ )
		for(iY = 2; iY < sim->ly; iY++)
		{
			if(sim->geometry[iX][iY] == FLUID)
			{
				int dX, dY;
				bool isNWC = false;
				
				for(dX = -1; dX <= 1; dX++ )
				{
					for(dY = -1; dY <= 1; dY++)
					{
						if(sim->geometry[iX+dX][iY+dY] == WALL || sim->geometry[iX+dX][iY+dY] == COAGWALL)
						{
							isNWC = true;
							//dY = 2;
							break;
						}
					}	
					if(isNWC)
						//dX = 2;
						break;
				}
				if(isNWC) {
					CellHistory* ch = createCellHistory();
					ch->x = iX; ch->y = iY;
					ch->rhoADP = defaultRhoADP;
					
					if(!arraylist_append(nwc, (void*)ch)) {
						printf("Memory allocation error!\n");
						exit(-1);
					}
					
					setGeometry(sim, iX, iY, NWFLUID);
				}
			}
		}
		
	return nwc;
}

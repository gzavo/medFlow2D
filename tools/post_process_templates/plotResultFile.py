# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 17:37:06 2015

@author: Gábor Závodszky

@description: This script processes a single result txt file.
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import os.path

# Read in the result files from here
workDir = 'C:/Data/temp/results/run1/'

# Value from simulation (e.g., the monitor file contains this),
# to calculate the physical time from the iteration number
dt = 3.576462e-5

# Set which file to plot from
fileName = 'full_0008945.txt'

# Load in the text file
def loadDataFile(fileName):

    fin = open(fileName)

    # Read in the first line with the geometry extents (resolution) of the output
    x,y = map(int, fin.readline().split())

    print " -> Loading result file ( %d x %d)..." % (x,y)
    
    # Select which values to process from the output file
    # Uncomment those necessary
    u=np.zeros((x,y))
    v=np.zeros((x,y))
    mag=np.zeros((x,y))
#    cX=np.zeros((x,y))
#    cY=np.zeros((x,y))
#    rho=np.zeros((x,y))
#    pres=np.zeros((x,y))
#    conc=np.zeros((x,y))
#    fx=np.zeros((x,y))
#    fy=np.zeros((x,y))

    
    for j in range(0, y):
        for i in range(0,x):
            # Read in the values from the text output file
            tu, tv, trho, tpres, tconc, tfx, tfy = map(float, fin.readline().split())

            # Uncomment the values you wish to store
            u[i,j]=tu
            v[i,j]=-tv
            mag[i,j]=np.sqrt(tu**2+tv**2)
#            pres[i,j]=tpres
#            rho[i,j]=trho
#            conc[i,j]=tconc
#            cX[i,j]=i
#            cY[i,j]=y-j
#            fx[i,j]=tfx
#            fy[i,j]=-tfy
        

    fin.close()

    # Return the desired results as a tuple
    return (x, y, u, v, mag)
    

# Change to the given directory
os.chdir(workDir)

# Call the loadDataFile function and convert the results to numpy arrays
x,y,u,v,mag = map(np.array, loadDataFile(fileName))

# Do some plotting with the result               
plt1 = plt.pcolor(mag.T, cmap=plt.get_cmap('jet'), vmin = 0, vmax = 0.11)
ax = plt.gca()
ax.set_aspect('equal')
ax.invert_yaxis()
plt.colorbar(plt1)

# Show the plot
plt.show()




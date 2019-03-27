# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 16:22:19 2015

@author: Gábor Závodszky

@description: This script processes an output folder populated with the txt output files.
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

# Save the computed values to this file in the result directory
outFileName = 'Output.txt'

# Save output plot to this file
outImageFile = 'Output.png'

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

    ### PUT YOUR CALCULATION WITH THE ARRAYS HERE

    # E.g., calculate the renal outflow in the abdominal geometry
    magLeft = 0
    for i in range(427,446+1):
        magLeft = magLeft + mag[14, i]
        
    magRight = 0
    for i in range(433, 454+1):
        magRight = magRight + mag[455, i]

    ### END OF CALCULATION 

    # Return the desired results as a tuple
    return (magLeft, magRight)



# Change to the given directory
os.chdir(workDir)

# Create some arrays to hold our calculations
# One value per time step
timeArray = []
leftMag = []
rightMag = []

# If the output file already exists, just read it back
if os.path.isfile(outFileName):
	with open(outFileName, "r") as text_file:
         for line in iter(text_file):
             t,l,r = map(float, line.split())
             timeArray.append(t)
             leftMag.append(l)
             rightMag.append(r)

# Otherwise create it, and calculate
else:    
    for fileName in glob.glob("*.txt"):
        s1 = fileName.find("_")+1
        s2 = fileName.find('.')

        if s1 < 1:      #if its a different kind of text file with no '_' skip it
            continue
    
        cTime = int(fileName[s1:s2])*dt
        print "Processing time:", cTime
    
        # Call the function to process a single text file
        lM,rM = loadDataFile(fileName)

        # Store the results
        timeArray.append(cTime)
        leftMag.append(lM)
        rightMag.append(rM)
	
    # Save the results to output.txt
    with open(outFileName, "w") as text_file:
        for j in range(0, len(timeArray)):
            text_file.write("%f %f %f\n" % (timeArray[j],leftMag[j], rightMag[j]))


# Do some plotting too with the results        
plt.plot(timeArray, leftMag)
plt.plot(timeArray, rightMag)

plt.legend(['left renal artery', 'right renal artery'], loc='upper left')

# Show the plot
# plt.show()

# Or rather save it as a png file
plt.savefig(outImageFile, bbox_inches='tight')
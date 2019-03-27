import numpy as np

fileName = "points.txt"

numPoints = 10000

fixIdx = 510.
fromIdx = 382.
toIdx = 386. 
#toIdx = 262.

rangeIdx = toIdx-fromIdx
step = rangeIdx / numPoints

with open(fileName, 'w') as outfile:
    for i in range(numPoints):
        yc = fixIdx
        xc = fromIdx + step*i
        outfile.write("%lf %lf\n" % (xc, yc))

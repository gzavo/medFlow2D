#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014, All Right Reserved, Gábor Závodszky, gabor@zavodszky.com
#
# This source is subject to the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# Please see the License.txt file for more information.
# All other rights reserved.
# You should have received a copy of the license along with this
# work in the License.txt file.
# If not, see <http://creativecommons.org/licenses/by-nc-nd/4.0/>.

import numpy as np
import matplotlib.pyplot as plt

x = []
y = []

coagProb = []
rhoADP = [] 
rhoPlateletAvg = []
shearMaxAvg = []

if __name__ == "__main__":
    import sys
    if(len(sys.argv) != 2):
        print "Argument needs to be a valid data file."
        sys.exit(-1)

    
    with open(sys.argv[1], "r") as f:
        for line in f:
            v1, v2, v3, v4, v5, v6 = [0,0,0,0,0,0]
            line_t = line.split()

            if line[0] == '#':
                if line[1] == '#':
                    #Not yet ready, needs more time in history queue
                    line_t.pop(0)
                    #!!!!!!!!!!!!!!!!!!
                    #If you do not wish to plot these values, uncomment the following 'continue'!
                    #!!!!!!!!!!!!!!!!!!

                    #continue                    
                else:
                    #It is a commented line
                    continue    

            try:
                v1, v2, v3, v4, v5, v6 = map(float, line_t)
            except ValueError:
                print "Line skipped (value error):", line
                continue

            x.append(int(v1))
            y.append(int(v2))
            coagProb.append(v3)
            rhoADP.append(v4)
            rhoPlateletAvg.append(v5)
            shearMaxAvg.append(v6)
    
    print len(x)

    plt.scatter(x,y, marker="s", c=coagProb, linewidths=0 )
    cax = plt.gca()
    cax.set_aspect('equal')
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.show()
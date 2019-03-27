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

SKIP = 10

u=None
v=None
mag=None
magF=None
cX=None
cY=None
rho=None
pres=None
conc=None
fx=None
fy=None

x=0
y=0

def loadDataFile(fileName):
    global u,v,mag,cX,cY,rho,pres,conc,fx,fy,magF,x,y

    fin = open(fileName)

    x,y = map(int, fin.readline().split())

    u=np.zeros((x,y))
    v=np.zeros((x,y))
    mag=np.zeros((x,y))
    magF=np.zeros((x,y))
    cX=np.zeros((x,y))
    cY=np.zeros((x,y))
    rho=np.zeros((x,y))
    pres=np.zeros((x,y))
    conc=np.zeros((x,y))
    fx=np.zeros((x,y))
    fy=np.zeros((x,y))

    
    for j in range(0, y):
        for i in range(0,x):
            tu, tv, trho, tpres, tconc, tfx, tfy = map(float, fin.readline().split())
            u[i,j]=tu
            v[i,j]=-tv
            mag[i,j]=np.sqrt(tu**2+tv**2)
            pres[i,j]=tpres
            rho[i,j]=trho
            conc[i,j]=tconc
            cX[i,j]=i
            cY[i,j]=y-j
            fx[i,j]=tfx
            fy[i,j]=-tfy 
            magF[i,j]=np.sqrt(tfx**2+tfy**2)

    fin.close()


if __name__ == "__main__":
    import sys
    if(len(sys.argv) != 3):
        print "Usage: python plotVector.py <datafile.txt> {0/1 - 0 for velocity, 1 for margination}"
        sys.exit(-1)

    loadDataFile(sys.argv[1])

    f1 = np.max(fx); f2 = np.max(fy)
    scaleF = (1.0 / max(f1,f2))
    
    v1 = np.max(u); v2 = np.max(v)
    scaleV = (1.0/max(v1,v2))

    plt.pcolor(cX,cY,mag)
    plt.colorbar()
    
    plt.contour(cX, cY, mag, [0,], linewidths=4, linestyles='solid', colors=('k',))

    ax = plt.gca()
    ax.set_aspect('equal')
    
    
    if(int(sys.argv[2])==1):
        #plt.quiver(cX[::SKIP, ::SKIP], cY[::SKIP, ::SKIP], fx[::SKIP, ::SKIP]*(1.0/magF[::SKIP, ::SKIP]), fy[::SKIP, ::SKIP]*(1.0/magF[::SKIP, ::SKIP]), width=0.001)
        plt.quiver(cX[::SKIP, ::SKIP], cY[::SKIP, ::SKIP], fx[::SKIP, ::SKIP]*scaleF, fy[::SKIP, ::SKIP]*scaleF, width=0.001)
        plt.title("Margination force [kg*m/s^2]")
    else:
        #plt.quiver(cX[::SKIP, ::SKIP], cY[::SKIP, ::SKIP], u[::SKIP, ::SKIP]*(1.0/mag[::SKIP, ::SKIP]), v[::SKIP, ::SKIP]*(1.0/mag[::SKIP, ::SKIP]), width=0.001)
        plt.quiver(cX[::SKIP, ::SKIP], cY[::SKIP, ::SKIP], u[::SKIP, ::SKIP], v[::SKIP, ::SKIP], width=0.001)
        plt.title("Velocity in [m/s]")

    plt.show()
  
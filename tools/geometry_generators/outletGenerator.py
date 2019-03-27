#!/usr/bin/python
# -*- coding: utf-8 -*-

from PIL import Image
import numpy as np
from tkFileDialog import askopenfilename
import Tkinter as tk

def Pixelreader(s1, wid, s2, he, bc):
    global o1, o2, in1, in2, tr1, tr2

    if bc == 0:
        o1 = []
        o2 = []
        in1 = []
        in2 = []
    elif bc == 1:
        tr1 = []
        tr2 = []

    for i in range(s1, wid, 1):
        for j in range(s2, he, 1):
            z = (j,i)
            r,g,b = img.getpixel((i,j))
            if ((r==0 and g==255 and b==0) and bc==0):
                o1 = np.append(o1, i)
                o2 = np.append(o2, j)
            elif ((r==255 and g==0 and b==0) and bc==0):
                in1 = np.append(in1, i)
                in2 = np.append(in2, j)
            elif (((r==0 or r==255) and g==255 and b==255) and bc==1):
                tr1 = np.append(tr1, i)
                tr2 = np.append(tr2, j)

#Opening the geometry file *.ppm

root = tk.Tk()
filename = askopenfilename()
img = Image.open(filename)
(width, height) = img.size
root.destroy()
print 'Opened file: %s, size: %dx%d \n' %(filename, width, height)

dinlet = int(raw_input("Enter the distance(in pixel) from the inlet:\n"))

outlet = []
inlet = []
o1 = []
o2 = []
in1 = []
in2 = []

#Reading green and red pixel coordinates

Pixelreader(0, width, 0, height, 0)

outlet = np.column_stack((o1,o2)) 
omin1 = int(outlet[0,0])
omin2 = int(outlet[0,1])
omax1 = int(outlet[-1,0]+1)
omax2 = int(outlet[-1,1]+1)
print '###########################################'
print 'Outlet: from (%d %d) to (%d %d) '% (omin1, omin2, omax1, omax2)

#with open("outlets.txt", 'w') as outfile:
#    outfile.write("%d %d %d %d" % (omin1, omin2, omax1, omax2))

inlet = np.column_stack((in1,in2)) 
imin1 = int(inlet[0,0])
imin2 = int(inlet[0,1])
imax1 = int(inlet[-1,0])
imax2 = int(inlet[-1,1])

print 'Inlet: from (%d %d) to (%d %d) '% (imin1, imin2, imax1, imax2)

trace1 = []
trace2 = []
if all(b == in1[0] for b in in1): # inlet is vertical
    if (in1[0] == 0): #left side
        print 'inlet: left side, vertical'
        trace1 = in1 + dinlet 
        trs1 = trace1[0]
        trwidth =  trs1 + 1
        if (in2[0] > 3 and in2[-1] < (height-3)):
            trs2 = in2[0] - 3
            trheight = in2[-1] + 3
        else:
            trs2 = in2[0]
            trheight = in2[-1] + 1

    else: #right side
        print 'inlet: right side, vertical'
        trace1 = in1 - dinlet
        trs1 = trace1[0]
        trwidth =  trs1+1
        if (in2[0] > 3 and in2[-1] < (height-3)):
            trs2 = in2[0]-3
            trheight = in2[-1] + 3
        else:
            trs2 = in2[0]
            trheight = in2[-1] + 1

elif all(bb == in2[0] for bb in in2): # inlet is horizontal
    if (in2[0] == 0): #top side
        print 'inlet: top side, horizontal'
        trace2 = in2 + dinlet 
        trs2 = trace2[0]
        trheight =  trs2+1
        if (in1[0] > 3 and (width - in1[-1]) > 3):
            trs1 = in1[0]-3
            trwidth = in1[-1] + 3
        else:
            trs1 = in1[0]
            trwidth = in1[-1] + 1
    else: # bottom side
        print 'inlet: bottom side, horizontal'
        trace2 = in2 - dinlet
        trs2 = trace2[0]
        trheight =  trs2 + 1
        if (in1[0] > 3 and (width - in1[-1]) > 3 ):
            trs1 = in1[0] - 3
            trwidth = in1[-1] + 3
        else:
            trs1 = in1[0]
            trwidth = in1[-1] + 1

trs1 = int(trs1)
trwidth = int(trwidth)
trs2 = int(trs2)
trheight = int(trheight)

#Reading the fluid coordinates at the start of the tracing

Pixelreader(trs1, trwidth, trs2, trheight , 1)

start = np.column_stack((tr1,tr2)) 
stmin1 = int(start[0,0])
stmin2 = int(start[0,1])
stmax1 = int(start[-1,0])
stmax2 = int(start[-1,1])

print 'Tracer points: from (%d %d) to (%d %d) '% (stmin1, stmin2, stmax1, stmax2)
print '###########################################'
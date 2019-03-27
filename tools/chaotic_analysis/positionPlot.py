#!/usr/bin/python
# -*- coding: utf-8 -*-

from tkFileDialog import askopenfilename
import numpy as np
import matplotlib.pyplot as plt
import Tkinter as tk
from PIL import Image
import sys

print '##################################################################'
print 'First open a geometry file *.ppm, then open a position file *.txt.'

root = tk.Tk()
root.withdraw()

filename1 = askopenfilename()
print 'Opened geometry: %s' %filename1
filename2 = askopenfilename()
print 'Opened pos file: %s' %filename2 
print '##################################################################'

print 'Choose one:\n 0 - load all points \n 1 - load without stuck points \n 2 - load only stuck points \n'
selectpoints = int(raw_input('Selected item: '))

if selectpoints == 1 or selectpoints == 2:
    print 'Open a res file *.txt'
    filename3 = askopenfilename()
    print 'Opened res file: %s' %filename3
    trtime = np.loadtxt(filename3, dtype = float, usecols=(0,))

print 'Working...'

im = Image.open(filename1)
(width, height) = im.size
im = im.transpose(Image.FLIP_TOP_BOTTOM)

xpos, ypos = np.loadtxt(filename2, dtype = float, usecols=(0,1), unpack = True)

if selectpoints == 0:
    newx = xpos
    newy = ypos
elif selectpoints == 1:
    k=0
    newx_temp = np.zeros_like(xpos)
    newy_temp = np.zeros_like(xpos)
    for i in range(len(trtime)):
        if trtime[i] != -1:
            newx_temp[k] = xpos[i]
            newy_temp[k] = ypos[i]
        k=k+1
    newx = newx_temp[0:k]
    newy = newy_temp[0:k]
elif selectpoints == 2:
    h=0
    newx_temp = np.zeros_like(xpos)
    newy_temp = np.zeros_like(xpos)
    for j in range(len(trtime)):
        if trtime[j] == -1:
            newx_temp[h] = xpos[j]
            newy_temp[h] = ypos[j]
        h=h+1
    newx = newx_temp[0:h]
    newy = newy_temp[0:h]
else:
    sys.exit('Error: the entered number is not 0, 1 or 2')

newy = height - newy
plt.imshow(im, origin='lower',extent=[-0.5, width-0.5, -0.5, height-0.5], interpolation='none') #extent=[0, width, 0, height]
plt.scatter(newx, newy, facecolor='r', marker='.', s=10.0, edgecolor='none')
plt.xlim(0, width)
plt.ylim(0, height)
plt.axes().set_aspect('equal')
plt.tick_params(which = 'both', direction='out')

print 'Done'
plt.show()
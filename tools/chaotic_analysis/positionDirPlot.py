#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from tkFileDialog import askopenfilename, askdirectory
import numpy as np
import matplotlib.pyplot as plt
import Tkinter as tk
from PIL import Image, ImageOps
import sys

print '##################################################################'
print 'First open a geometry file *.ppm, then open a directory that contains position files *.txt.'

root = tk.Tk()
root.withdraw()

filename1 = askopenfilename()
filename1 = os.path.normpath(filename1)
print 'Opened geometry: %s' %filename1

filename2 = askdirectory()
filename2 = os.path.normpath(filename2)
filename2 = os.path.join(filename2,'')
pfilelist = os.listdir(filename2)
posfiles = []
for pfile in pfilelist:
    if pfile.endswith(".txt"):
        if pfile.startswith("res_pos"):
            posfiles.append(pfile)

print 'Opened pos directory: %s' %filename2 
print '##################################################################'

print 'Choose one:\n 0 - load all points \n 1 - load without stuck points \n 2 - load only stuck points \n'
selectpoints = int(raw_input('Selected item: '))

if selectpoints == 1 or selectpoints == 2:
    print 'Open res files directory'
    filename3 = askdirectory()
    filename3 = os.path.normpath(filename3)
    filename3 = os.path.join(filename3,'')
    rfilelist = os.listdir(filename3)
    resfiles = []
    for rfile in rfilelist:
        if rfile.endswith(".txt"):
            if rfile.startswith("res_res"):
                resfiles.append(rfile)
    print 'Opened res directory: %s' %filename3

print 'Working...'

im = Image.open(filename1)
(width, height) = im.size
im = im.transpose(Image.ROTATE_180)
im = im.transpose(Image.FLIP_LEFT_RIGHT)

out = filename2 +'imgpos'
out = os.path.normpath(out)
out = os.path.join(out,'')

print 'Saving images to: %s' %out
if not os.path.exists(out):
    os.makedirs(out)

for j in posfiles:
    x, y = [0, 0]
    resi = 0
    xpos, ypos = np.loadtxt(filename2+j, dtype = float, usecols=(0,1), unpack = True)
    if selectpoints == 0:
        newx = xpos
        newy = ypos
    elif selectpoints == 1:
        k=0
        trtime = np.loadtxt(filename3+resfiles[resi], dtype = float, usecols=(0,))
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
        trtime = np.loadtxt(filename3+resfiles[resi], dtype = float, usecols=(0,))
        newx_temp = np.zeros_like(xpos)
        newy_temp = np.zeros_like(xpos)
        for n in range(len(trtime)):
            if trtime[n] == -1:
                newx_temp[h] = xpos[n]
                newy_temp[h] = ypos[n]
            h=h+1
        newx = newx_temp[0:h]
        newy = newy_temp[0:h]
    else:
        sys.exit('Error: the entered number is not 0, 1 or 2')

    newy = height-newy
    plt.imshow(im, origin='lower', extent=[-0.5, width-0.5, -0.5, height-0.5], interpolation="none")
    plt.scatter(newx, newy, facecolor='r', marker=".", s=0.5, edgecolor='none')
    plt.xlim(0, width)
    plt.ylim(0, height)
    plt.axes().set_aspect('equal')
    plt.tick_params(which = 'both', direction='out')
    imgname = '%s.png' %j[:-5]
    plt.savefig(out+imgname, dpi=1000, bbox_inches='tight')
    plt.close()
    resi = resi+1
    print 'Saving: %s' %imgname
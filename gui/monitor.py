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
from __future__ import division
from __future__ import print_function

#import tkFileDialog
#from Tkinter import *
try:
    from Tkinter import *
except ImportError:
    from tkinter import *

try:
    import ttk
    from ttk import Notebook
    py3 = 0
except ImportError:
    import tkinter.ttk as ttk
    from tkinter.ttk import Notebook
    py3 = 1

if py3:
    from tkinter import filedialog as tkFileDialog
else:
    import tkFileDialog

from PIL import Image, ImageTk #import Image, ImageTk
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
import sys
import os
import math

import zipfile

baseDir = '..'

image_name1 = '/monitor/currentVel.png'
image_name2 = '/monitor/currentRho.png'
image_name3 = '/monitor/currentCoupledRho.png'
text_name   = '/monitor/currentInfo.txt'
text_vel_name='currentVel.txt'
text_vel_archive_name = '/monitor/'+text_vel_name+'.zip'

updInterval = 2500

inDir = None

scaleF = 0

class AutoScrollbar(Scrollbar):
    # a scrollbar that hides itself if it's not needed.  only
    # works if you use the grid geometry manager.
    # autoScrollbar code from http://effbot.org/zone/tkinter-autoscrollbar.htm
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            # grid_remove is currently missing from Tkinter!
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)
    def pack(self, **kw):
        raise TclError("cannot use pack with this widget")
    def place(self, **kw):
        raise TclError ("cannot use place with this widget")

def loadText(fileName):
    with open (fileName, "r") as myfile:
        return myfile.read()
 
def update_image():
    global tkimg1, image_name1, tkimg2, image_name2, tkimg3, image_name3, updInterval
    
    try:
        tkimg1 = ImageTk.PhotoImage(Image.open(baseDir+image_name1))
        label.config( image = tkimg1)
        tkimg2 = ImageTk.PhotoImage(Image.open(baseDir+image_name2))
        label2.config( image = tkimg2)
        tkimg3 = ImageTk.PhotoImage(Image.open(baseDir+image_name3))
        label3.config( image = tkimg3)
        label4.config( text=loadText(baseDir+text_name))
        label.after(updInterval, update_image)
    #except Exception:
    except (IOError, SyntaxError):
        label.after(updInterval, update_image)
        #print "Image file is not ready, update postponed."
        

def choose_dir():
    global baseDir, inDir
    baseDir = tkFileDialog.askdirectory()
    inDir.config(text=baseDir)
    print("Base directory change: " + baseDir)

def plotScaled():
    import matplotlib.pyplot as plt
    import numpy as np

    zip_fin = zipfile.ZipFile(baseDir+text_vel_archive_name)
    fin = zip_fin.open(text_vel_name)

    x,y = map(int, fin.readline().split())

    u=np.zeros((x,y))
    v=np.zeros((x,y))
    mag=np.zeros((x,y))
    cX=np.zeros((x,y))
    cY=np.zeros((x,y))
    rho=np.zeros((x,y))
    pres=np.zeros((x,y))
    conc=np.zeros((x,y))
    fx=np.zeros((x,y))
    fy=np.zeros((x,y))
    
    try:

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
    except ValueError:
        fin.close()
        return

    fin.close()
    zip_fin.close()

    f1 = np.max(fx); f2 = np.max(fy)
    forcemax = max(f1,f2)
    if forcemax != 0:
        scaleF = (1.0 / forcemax )*0.025
    else:
        scaleF = 1

    v1 = np.max(u); v2 = np.max(v)
    velmax = max(v1,v2)
    if velmax !=0:
        scaleV = (1.0/velmax )*0.025
    else:
        scaleV = 1

    plt.figure(1)
    #plt.subplot(211)
    if CheckVar1.get():
        mgmin = mag[mag !=0].min()
        mgmax = mag.max()
        lgmgmin = math.log10(mgmin)
        lgmgmax = math.log10(mgmax)
        plt.pcolor(cX,cY, mag, norm=LogNorm(vmin=mgmin, vmax=mgmax ))
        ax = plt.gca()
        ax.set_aspect('equal')
        lvls = np.logspace(lgmgmin, lgmgmax, 10)
        l_f = LogFormatter(10, labelOnlyBase=False)
        plt.colorbar(ticks=lvls, format=l_f)
        plt.quiver(cX, cY, u*scaleV, v*scaleV, scale=0.4, width=0.001)
        plt.title("Log scaled velocity in [m/s]")
    else:
        plt.pcolor(cX,cY,mag)
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.colorbar()
        plt.quiver(cX, cY, u*scaleV, v*scaleV, scale=0.4, width=0.001)
        plt.title("Velocity in [m/s]")
    
    plt.figure(2)
    #plt.subplot(212)
    plt.pcolor(cX,cY,pres)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.colorbar()
    plt.title("Relative pressure in [Pa]")
    
    if not np.all(conc == 0):
        plt.figure(3)
        #plt.subplot(211)
        plt.pcolor(cX,cY,conc)
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.colorbar()
        if not (np.all(fx == 0) and np.all(fy == 0)):
            plt.quiver(cX, cY, fx*scaleF, fy*scaleF, scale=0.4, width=0.001)
            plt.title("Concentration and acting forces [kg*m/s^2]")
        else:
            plt.title("Fluid age [s]")
    
    plt.show()


if __name__ == "__main__":
    
    if(len(sys.argv) == 2):
        if(os.path.isfile(sys.argv[1])):
            if py3:
                import configparser
                config = configparser.ConfigParser(inline_comment_prefixes=';')
            else:
                import ConfigParser
                config = ConfigParser.ConfigParser()
        
            
            config.read(sys.argv[1])
            baseDir=config.get('output','workingDir')
            
        elif(os.path.isdir(sys.argv[1])):
            baseDir=sys.argv[1]
            
        else:
            print("Argument needs to be either a valid directory or an input ini file.")
            sys.exit(-1)
                 
    else:
        print ("Usage:", sys.argv[0], " <setup.ini file or output directory>" )
        sys.exit(-1)
        
        
    w = Tk()
    w.wm_title("medFlow2D monitor")
    #--------------------------------------------------------------------------------
    vscrollbar = AutoScrollbar(w)
    vscrollbar.grid(row=0, column=1, sticky=N+S)
    hscrollbar = AutoScrollbar(w, orient=HORIZONTAL)
    hscrollbar.grid(row=1, column=0, sticky=E+W)
    canvas = Canvas(w, width=1024, height=400,
                    yscrollcommand=vscrollbar.set,
                    xscrollcommand=hscrollbar.set,
                    highlightthickness=0)
    canvas.grid(row=0, column=0, sticky=N+S+E+W)
 
    vscrollbar.config(command=canvas.yview)
    hscrollbar.config(command=canvas.xview)

    # make the canvas expandable
    w.grid_rowconfigure(0, weight=1)
    w.grid_columnconfigure(0, weight=1)

    #
    # create canvas contents

    frame = Frame(canvas)

    frame.rowconfigure(1, weight=1)
    frame.columnconfigure(1, weight=1)
    
    #--------------------------------------------------------------------------------
    
    b1 = Button(frame, text="Plot detailed fields",command=plotScaled)
    b1.grid(row=0, column=1, sticky='W')

    CheckVar1 = IntVar()
    b2 = Checkbutton(frame, text = "Log scaled velocity", variable = CheckVar1, onvalue = 1, offvalue = 0)
    b2.grid(row=0, column=0, sticky='E')
  
    b3 = Button(frame, text="Switch base dir",command=choose_dir)
    b3.grid(row=0, column=3, sticky='E')

    inDir = Label(frame, text=baseDir)
    inDir.grid(row=0, column=4, sticky='W')

    #Notebook widget
    note = Notebook(frame)
    note.grid(row=1, column=0, columnspan=4, sticky='W')
    
    #First image
    im = Image.open(baseDir+image_name1)
    tkimg1 = ImageTk.PhotoImage(im)  
    label =  Label(frame, image=tkimg1)
    label.grid(row=1, column=0, columnspan=4 )
    note.add(label, text='Velocity')

    #Second image
    im2 = Image.open(baseDir+image_name2)
    tkimg2 = ImageTk.PhotoImage(im2)
    label2 =  Label(frame, image=tkimg2)
    label2.grid(row=1, column=0, columnspan=4  )
    note.add(label2, text='Density')

    #Third image
    im3 = Image.open(baseDir+image_name3)
    tkimg3 = ImageTk.PhotoImage(im3)
    label3 =  Label(frame, image=tkimg3)
    label3.grid(row=1, column=0, columnspan=4  )
    note.add(label3, text='Coupled field')

    #Textual data
    label4 =  Label(frame, text=loadText(baseDir+text_name), justify='left', padx=5)
    label4.grid(row=0, column=4, rowspan=2 )
    
    #-------------------------------------------------------------------------------- 
    
    canvas.create_window(0, 0, anchor=NW, window=frame)

    frame.update_idletasks()

    canvas.config(scrollregion=canvas.bbox("all"))
    
    w.after(updInterval, update_image)
    
    w.wm_iconbitmap('@gui/main.xbm')
    w.mainloop()
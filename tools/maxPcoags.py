#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
import numpy as np

import matplotlib.pyplot as plt
import Tkinter, tkFileDialog, Tkconstants
from Tkinter import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

baseDir = 'C:/'
coagfiles=[]  
iteration=[]
Pcoags=[]
nerror=0
coagmax=[]
intiteration=[]

def choose_dir():
    global baseDir, L1, Text1, Text2
    baseDir = tkFileDialog.askdirectory()
    baseDir= os.path.normpath(baseDir)
    baseDir = os.path.join(baseDir,'')
    Text1 = "Base directory:  %s" % baseDir
    Text2 = "  ** Press Load!**"
    print  "%s   %s" % (Text1, Text2)
    L1.config(text="%s   %s" % (Text1, Text2))

def _quit():
    global root
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent

                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
def SaveData():
    sfileName = tkFileDialog.asksaveasfile(mode='w', defaultextension='txt')

    if(sfileName):        
        for i in range(len(iteration)):
            sfileName.write(str(iteration[i])+"  "+str(coagmax[i])+"\n")
        sfileName.close()

def on_key_event(event):
    global canvas, toolbar
    print('you pressed %s'%event.key)
    key_press_handler(event, canvas, toolbar)

def ClearPlot():
    global f, L1, Text1, Text2
    ax.cla()
    ax.set_yscale('log')
    ax.set_xlabel('Number of iteration')
    ax.set_ylabel('Pcoag')
    canvas.draw()
    Text2 = "  ** Plot cleared!**"
    print  "%s   %s" % (Text1, Text2)
    L1.config(text="%s   %s" % (Text1, Text2))
    
def plotPcoag():
    global f, L1, Text1, Text2
    graph1 = f.add_subplot(111)
    graph1.plot(iteration,coagmax)
    canvas.draw()
    Text2 = "  ** Plot added!**"
    print  "%s   %s" % (Text1, Text2)
    L1.config(text="%s   %s" % (Text1, Text2))


def load_data():
    global coagfiles, iteration, Pcoags, nerror, coagmax, intiteration, L1, Text1, Text2
    dirs = os.listdir(baseDir)
    coagfiles=[]
    for filename in dirs:
        if filename.endswith(".txt"):
            if filename.startswith("coag"):
                coagfiles.append(filename)
    nerror=0
    iteration=[]
    coagmax=[]
    intiteration=[]
    for coagfname in coagfiles:
        with open(baseDir+coagfname, "r") as f:
            line=0
            Pcoags=[]
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
                    nerror+=1
                    continue
            
                Pcoags.append(v3)
                niteration=str(coagfname[4:-4])
                i=0
                while (niteration[i]==0):
                    i+=i
                szam=''
                for i in xrange(1,len(niteration)):
                    szam=szam+niteration[i]

            intiteration=int(szam)
            iteration.append(intiteration)
            coagmax.append(max(Pcoags))

        f.close()     
    Text2='Loaded, Value error: %d,  **Press Plot!**' %nerror
    print  "%s   %s" % (Text1, Text2)
    L1.config(text="%s   %s" % (Text1, Text2))
    

if __name__=='__main__':
    root = Tk()

    root.configure(background='grey')
    f = Figure(figsize=(10,5), dpi=100)
    t=iteration[:]
    s=coagmax[:]

    ax = f.add_subplot(1,1,1) 
    ax.plot(t,s)
    
    ax.set_yscale('log')
    ax.set_xlabel('Number of iteration')
    ax.set_ylabel('Pcoag')
    
    
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=1)
    toolbar=NavigationToolbar2TkAgg(canvas,root)
    toolbar.pack(side=BOTTOM)
    toolbar.grid(row=2, column=0, sticky="ew")
    canvas.mpl_connect('key_press_event', on_key_event)

    toolbar = Frame(root)
    
    Text1 = "Base directory:  %s" % baseDir
    Text2 = "  ** Choose a dir!**"

    b1 = Button(root, text="Base dir",command=choose_dir)
    b1.grid(row=0, column=0,sticky=W, padx=5)

    b2=Button(root, text='Load', command=load_data)
    b2.grid(row=0, column=0, sticky=W, padx=60)

    b3=Button(root, text='Plot', command=plotPcoag)
    b3.grid(row=0, column=0, sticky=W, padx=100)

    b4=Button(root, text='Save', command=SaveData)
    b4.grid(row=0, column=0, sticky=W, padx=135)

    b6=Button(root, text='Clear', command=ClearPlot)
    b6.grid(row=0, column=0, sticky=W, padx=175)

    L1=Label(root, bg="white", text = "%s   %s" % (Text1,Text2))
    L1.grid(row=0, column=0, sticky=W, padx=220, columnspan=5)

    b5 = Button(master=root, text='Quit', command=_quit)
    b5.grid(row=0, column=5)

    toolbar.grid(row=3, column=0, sticky="ew") 

    canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew", columnspan=6)
    root.grid_rowconfigure(0, weight=0)
    root.grid_rowconfigure(1, weight=1)
    root.grid_columnconfigure(1, weight=1)    
    
    root.mainloop()
#! /usr/bin/env python
#
# GUI module generated by PAGE version 4.5
# In conjunction with Tcl version 8.6
#    Jun 22, 2015 02:46:40 PM
from __future__ import division
from __future__ import print_function
import sys

try:
    from Tkinter import *
except ImportError:
    from tkinter import *

try:
    import ttk
    py3 = 0
except ImportError:
    import tkinter.ttk as ttk
    py3 = 1

if py3:
    from tkinter import filedialog as tkFileDialog
else:
    import tkFileDialog

import solver_control_support

import os
import subprocess

from sys import platform as _platform
if _platform == "linux" or _platform == "linux2":
    binaryFileName = "/../bin/medFlow2D"
elif _platform == "darwin":
    binaryFileName = "/../bin/medFlow2D_mac"
elif _platform == "win32":
    binaryFileName = "/../bin/medFlow2D.exe"

monitorFileName = "/../gui/monitor.py"
editorFileName = "/../gui/editor.py"
plotFileName = "/../gui/plotData.py"

setupFileName = "Choose a setup file!"
postDirectory = "."

solverProc = None

numFiles = 0

strDt = None

def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root
    root = Tk()
    root.title('medFlow2D solver control')
    geom = "454x278+647+253"
    root.geometry(geom)
    w = medFlow2D_solver_control (root)
    solver_control_support.init(root, w)
    root.wm_iconbitmap('@gui/main.xbm')
    root.mainloop()

w = None
def create_medFlow2D_solver_control(root, param=None):
    '''Starting point when module is imported by another program.'''
    global w, w_win, rt
    rt = root
    w = Toplevel (root)
    w.title('medFlow2D solver control')
    geom = "454x278+647+253"
    w.geometry(geom)
    w_win = medFlow2D_solver_control (w)
    solver_control_support.init(w, w_win, param)
    return w_win

def destroy_medFlow2D_solver_control():
    global w
    w = None
    root.destroy()

def setupNext(notebook):
    global solverProc

    if _platform == "linux" or _platform == "linux2" or _platform == "darwin":
        solverProc = subprocess.Popen([os.path.dirname(os.path.realpath(__file__)) + binaryFileName + " " + setupFileName], shell=True)
    elif _platform == "win32":
        solverProc = subprocess.Popen([os.path.dirname(os.path.realpath(__file__)) + binaryFileName, setupFileName], shell=True)

    notebook.tab(1, state='normal')
    notebook.tab(0, state='disabled')

def saveSim():
    cmdDirectory = os.path.dirname(setupFileName)
    cmdFile = cmdDirectory + "/.command"

    outFile = open(cmdFile, 'w')
    outFile.write("SAVEA")
    outFile.close()


def stopSim(notebook):
    cmdDirectory = os.path.dirname(setupFileName)
    cmdFile = cmdDirectory + "/.command"

    try:
        outFile = open(cmdFile, 'w')
        outFile.write("EXIT")
        outFile.close()
    except IOError as e:
        print(str(e))

    notebook.tab(0, state='normal')
    notebook.tab(1, state='disabled')
    notebook.select(0)


def openSetup(label):
    global setupFileName
    options = {}
    options['defaultextension'] = '.ini'
    options['filetypes'] = [('ini files', '.ini')]
    options['initialdir'] = os.path.dirname(os.path.realpath(__file__)) + '../data'
    options['initialfile'] = 'setup.ini'
    options['parent'] = root
    options['title'] = 'Open simulation setup'
    newSetupFile = tkFileDialog.askopenfilename(**options)
    if(os.path.isfile(newSetupFile)):
        setupFileName = newSetupFile
        label.configure(text=setupFileName)

def startMonitor():
    #cmdMonitor = os.path.dirname(os.path.realpath(__file__)) + monitorFileName + " " + setupFileName
    if _platform == "linux" or _platform == "linux2" or _platform == "darwin":
        subprocess.Popen(["python", os.path.dirname(os.path.realpath(__file__)) + monitorFileName, setupFileName], shell=False)
    elif _platform == "win32":
        subprocess.Popen(["python", os.path.dirname(os.path.realpath(__file__)) + monitorFileName, setupFileName], shell=False)

def openEdit():
    #cmdEditor = os.path.dirname(os.path.realpath(__file__)) + editorFileName + " " + setupFileName
    if _platform == "linux" or _platform == "linux2" or _platform == "darwin":
        subprocess.Popen(["python", os.path.dirname(os.path.realpath(__file__)) + editorFileName, setupFileName], shell=False)
    elif _platform == "win32":
        subprocess.Popen(["python", os.path.dirname(os.path.realpath(__file__)) + editorFileName, setupFileName], shell=False)

def openPostDir(label, listbox):
    global postDirectory, numFiles
    options = {}

    postDir = tkFileDialog.askdirectory()

    if(os.path.isdir(postDir)):
        postDirectory = postDir
        label.configure(text=postDirectory)

        for i in range(0, numFiles):
            listbox.delete(0, 0)

        numFiles = 0
        for afile in os.listdir(postDir):
            if afile.endswith(".txt.zip"):
                listbox.insert(END, afile)
                numFiles = numFiles + 1

def reloadPostDir(listbox):
    global postDirectory, numFiles

    for i in range(0, numFiles):
        listbox.delete(0, 0)

    numFiles = 0
    for afile in os.listdir(postDirectory):
        if afile.endswith(".txt.zip"):
            listbox.insert(END, afile)
            numFiles = numFiles + 1

def openPlot(listbox):
    selected = listbox.curselection()
    if(not selected):
        return
    idx = selected[0]
    toPlotFile = postDirectory + '/' + listbox.get(idx)

    if _platform == "linux" or _platform == "linux2" or _platform == "darwin":
        subprocess.Popen(["python", os.path.dirname(os.path.realpath(__file__)) + plotFileName, toPlotFile], shell=False)
    elif _platform == "win32":
        subprocess.Popen(["python", os.path.dirname(os.path.realpath(__file__)) + plotFileName, toPlotFile], shell=False)

def calcTime(label, listbox, evt):
    selected = listbox.curselection()
    if(not selected):
        return
    idx = selected[0]
    strName = listbox.get(idx)

    #print float(strDt.get()), float(strName[5:12])
    label.configure( text="t = " + str(round(float(strName[5:12]) * float(strDt.get()), 5)) + " s" )


class medFlow2D_solver_control:
    def __init__(self, master=None):
        global strDt

        _bgcolor = '#d9d9d9'  # X11 color: 'gray85'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9' # X11 color: 'gray85'
        _ana1color = '#d9d9d9' # X11 color: 'gray85'
        _ana2color = '#d9d9d9' # X11 color: 'gray85'
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.',background=_bgcolor)
        self.style.configure('.',foreground=_fgcolor)
        self.style.configure('.',font="TkDefaultFont")
        self.style.map('.',background=
            [('selected', _compcolor), ('active',_ana2color)])
        master.configure(background="#d9d9d9")
        master.configure(highlightbackground="#d9d9d9")
        master.configure(highlightcolor="black")


        self.style.configure('TNotebook.Tab', background=_bgcolor)
        self.style.configure('TNotebook.Tab', foreground=_fgcolor)
        self.style.map('TNotebook.Tab', background=
            [('selected', _compcolor), ('active',_ana2color)])
        self.TNotebook1 = ttk.Notebook(master)
        self.TNotebook1.place(relx=0.02, rely=0.04, relheight=0.92
                , relwidth=0.96)
        self.TNotebook1.configure(width=434)
        self.TNotebook1.configure(takefocus="")
        self.TNotebook1_pg0 = ttk.Frame(self.TNotebook1)
        self.TNotebook1.add(self.TNotebook1_pg0, padding=3)
        self.TNotebook1.tab(0, text="Setup",underline="-1",)
        self.TNotebook1_pg1 = ttk.Frame(self.TNotebook1)
        self.TNotebook1.add(self.TNotebook1_pg1, padding=3)
        self.TNotebook1.tab(1, text="Run",underline="-1", state='disabled')
        self.TNotebook1_pg2 = ttk.Frame(self.TNotebook1)
        self.TNotebook1.add(self.TNotebook1_pg2, padding=3)
        self.TNotebook1.tab(2, text="Post",underline="-1",)

        self.Label1 = Label(self.TNotebook1_pg0)
        self.Label1.place(relx=0.02, rely=0.09, height=21, width=74)
        self.Label1.configure(activebackground="#f9f9f9")
        self.Label1.configure(activeforeground="black")
        self.Label1.configure(background=_bgcolor)
        self.Label1.configure(disabledforeground="#a3a3a3")
        self.Label1.configure(foreground="#000000")
        self.Label1.configure(highlightbackground="#d9d9d9")
        self.Label1.configure(highlightcolor="black")
        self.Label1.configure(text='''Setup ini file:''')

        self.Label2 = Label(self.TNotebook1_pg0)
        self.Label2.place(relx=0.02, rely=0.17, height=21, width=344)
        self.Label2.configure(activebackground="#f9f9f9")
        self.Label2.configure(activeforeground="black")
        self.Label2.configure(background=_bgcolor)
        self.Label2.configure(disabledforeground="#a3a3a3")
        self.Label2.configure(foreground="#000000")
        self.Label2.configure(highlightbackground="#d9d9d9")
        self.Label2.configure(highlightcolor="black")
        self.Label2.configure(justify=LEFT)
        self.Label2.configure(relief=SUNKEN)
        self.Label2.configure(text=setupFileName)

        self.Button1 = Button(self.TNotebook1_pg0)
        self.Button1.place(relx=0.84, rely=0.17, height=24, width=57)
        self.Button1.configure(activebackground="#d9d9d9")
        self.Button1.configure(activeforeground="#000000")
        self.Button1.configure(background=_bgcolor)
        self.Button1.configure(disabledforeground="#a3a3a3")
        self.Button1.configure(foreground="#000000")
        self.Button1.configure(highlightbackground="#d9d9d9")
        self.Button1.configure(highlightcolor="black")
        self.Button1.configure(pady="0")
        self.Button1.configure(text='''Open''')
        self.Button1.configure(command=lambda: openSetup(self.Label2))

        self.Button2 = Button(self.TNotebook1_pg0)
        self.Button2.place(relx=0.02, rely=0.39, height=24, width=407)
        self.Button2.configure(activebackground="#d9d9d9")
        self.Button2.configure(activeforeground="#000000")
        self.Button2.configure(background=_bgcolor)
        self.Button2.configure(disabledforeground="#a3a3a3")
        self.Button2.configure(foreground="#000000")
        self.Button2.configure(highlightbackground="#d9d9d9")
        self.Button2.configure(highlightcolor="black")
        self.Button2.configure(pady="0")
        self.Button2.configure(text='''Edit''')
        self.Button2.configure(command=openEdit)

        self.Button3 = Button(self.TNotebook1_pg0)
        self.Button3.place(relx=0.02, rely=0.57, height=24, width=407)
        self.Button3.configure(activebackground="#d9d9d9")
        self.Button3.configure(activeforeground="#000000")
        self.Button3.configure(background=_bgcolor)
        self.Button3.configure(disabledforeground="#a3a3a3")
        self.Button3.configure(foreground="#000000")
        self.Button3.configure(highlightbackground="#d9d9d9")
        self.Button3.configure(highlightcolor="black")
        self.Button3.configure(pady="0")
        self.Button3.configure(text='''Run simulation''')
        self.Button3.configure(command= lambda: setupNext(self.TNotebook1))

        self.Button4 = Button(self.TNotebook1_pg0)
        self.Button4.place(relx=0.02, rely=0.83, height=24, width=407)
        self.Button4.configure(activebackground="#d9d9d9")
        self.Button4.configure(activeforeground="#000000")
        self.Button4.configure(background=_bgcolor)
        self.Button4.configure(disabledforeground="#a3a3a3")
        self.Button4.configure(foreground="#000000")
        self.Button4.configure(highlightbackground="#d9d9d9")
        self.Button4.configure(highlightcolor="black")
        self.Button4.configure(pady="0")
        self.Button4.configure(text='''Exit''')
        self.Button4.configure(command=destroy_medFlow2D_solver_control)

        self.Label3 = Label(self.TNotebook1_pg1)
        self.Label3.place(relx=0.05, rely=0.09, height=21, width=54)
        self.Label3.configure(activebackground="#f9f9f9")
        self.Label3.configure(activeforeground="black")
        self.Label3.configure(background=_bgcolor)
        self.Label3.configure(disabledforeground="#a3a3a3")
        self.Label3.configure(foreground="#000000")
        self.Label3.configure(highlightbackground="#d9d9d9")
        self.Label3.configure(highlightcolor="black")
        self.Label3.configure(text='''Status:''')

        self.Button5 = Button(self.TNotebook1_pg1)
        self.Button5.place(relx=0.07, rely=0.26, height=24, width=377)
        self.Button5.configure(activebackground="#d9d9d9")
        self.Button5.configure(activeforeground="#000000")
        self.Button5.configure(background=_bgcolor)
        self.Button5.configure(disabledforeground="#a3a3a3")
        self.Button5.configure(foreground="#000000")
        self.Button5.configure(highlightbackground="#d9d9d9")
        self.Button5.configure(highlightcolor="black")
        self.Button5.configure(pady="0")
        self.Button5.configure(text='''Start a new monitor''')
        self.Button5.configure(command=startMonitor)

        self.Button6 = Button(self.TNotebook1_pg1)
        self.Button6.place(relx=0.07, rely=0.83, height=24, width=374)
        self.Button6.configure(activebackground="#d9d9d9")
        self.Button6.configure(activeforeground="#000000")
        self.Button6.configure(background=_bgcolor)
        self.Button6.configure(disabledforeground="#a3a3a3")
        self.Button6.configure(foreground="#000000")
        self.Button6.configure(highlightbackground="#d9d9d9")
        self.Button6.configure(highlightcolor="black")
        self.Button6.configure(pady="0")
        self.Button6.configure(text='''Stop simulation and go back to setup''')
        self.Button6.configure(command= lambda: stopSim(self.TNotebook1))

        self.Label5 = Label(self.TNotebook1_pg1)
        self.Label5.place(relx=0.21, rely=0.09, height=21, width=314)
        self.Label5.configure(activebackground="#f9f9f9")
        self.Label5.configure(activeforeground="black")
        self.Label5.configure(background=_bgcolor)
        self.Label5.configure(disabledforeground="#a3a3a3")
        self.Label5.configure(foreground="#000000")
        self.Label5.configure(highlightbackground="#d9d9d9")
        self.Label5.configure(highlightcolor="black")
        self.Label5.configure(relief=RIDGE)
        self.Label5.configure(text='''Simulation started.''')

        self.Button7 = Button(self.TNotebook1_pg1)
        self.Button7.place(relx=0.09, rely=0.48, height=24, width=167)
        self.Button7.configure(activebackground="#d9d9d9")
        self.Button7.configure(activeforeground="#000000")
        self.Button7.configure(background=_bgcolor)
        self.Button7.configure(disabledforeground="#a3a3a3")
        self.Button7.configure(foreground="#000000")
        self.Button7.configure(highlightbackground="#d9d9d9")
        self.Button7.configure(highlightcolor="black")
        self.Button7.configure(pady="0")
        self.Button7.configure(text='''Create checkpoint now''')
        self.Button7.configure(width=167)

        self.Button10 = Button(self.TNotebook1_pg1)
        self.Button10.place(relx=0.09, rely=0.61, height=24, width=347)
        self.Button10.configure(activebackground="#d9d9d9")
        self.Button10.configure(activeforeground="#000000")
        self.Button10.configure(background=_bgcolor)
        self.Button10.configure(disabledforeground="#a3a3a3")
        self.Button10.configure(foreground="#000000")
        self.Button10.configure(highlightbackground="#d9d9d9")
        self.Button10.configure(highlightcolor="black")
        self.Button10.configure(pady="0")
        self.Button10.configure(text='''Save all data now''')
        self.Button10.configure(width=347)
        self.Button10.configure(command=saveSim)

        self.Button11 = Button(self.TNotebook1_pg1)
        self.Button11.place(relx=0.51, rely=0.48, height=24, width=167)
        self.Button11.configure(activebackground="#d9d9d9")
        self.Button11.configure(activeforeground="#000000")
        self.Button11.configure(background=_bgcolor)
        self.Button11.configure(disabledforeground="#a3a3a3")
        self.Button11.configure(foreground="#000000")
        self.Button11.configure(highlightbackground="#d9d9d9")
        self.Button11.configure(highlightcolor="black")
        self.Button11.configure(pady="0")
        self.Button11.configure(text='''Save images now''')
        self.Button11.configure(width=167)

        self.Scrolledlistbox1 = ScrolledListBox(self.TNotebook1_pg2)
        self.Scrolledlistbox1.place(relx=0.02, rely=0.17, relheight=0.78
                , relwidth=0.72)
        self.Scrolledlistbox1.configure(background="white")
        self.Scrolledlistbox1.configure(disabledforeground="#a3a3a3")
        self.Scrolledlistbox1.configure(font="TkFixedFont")
        self.Scrolledlistbox1.configure(foreground="black")
        self.Scrolledlistbox1.configure(highlightbackground="#d9d9d9")
        self.Scrolledlistbox1.configure(highlightcolor="#d9d9d9")
        self.Scrolledlistbox1.configure(selectbackground="#c4c4c4")
        self.Scrolledlistbox1.configure(selectforeground="black")
        self.Scrolledlistbox1.configure(width=30)
        self.Scrolledlistbox1.configure(selectmode=SINGLE)
        #self.Scrolledlistbox1.configure(command=lambda: calcTime(self.Label7, self.Scrolledlistbox1))
        self.Scrolledlistbox1.bind('<<ListboxSelect>>', lambda e: calcTime(self.Label7, self.Scrolledlistbox1, e))

        self.Button8 = Button(self.TNotebook1_pg2)
        self.Button8.place(relx=0.77, rely=0.04, height=24, width=87)
        self.Button8.configure(activebackground="#d9d9d9")
        self.Button8.configure(activeforeground="#000000")
        self.Button8.configure(background=_bgcolor)
        self.Button8.configure(disabledforeground="#a3a3a3")
        self.Button8.configure(foreground="#000000")
        self.Button8.configure(highlightbackground="#d9d9d9")
        self.Button8.configure(highlightcolor="black")
        self.Button8.configure(pady="0")
        self.Button8.configure(text='''Open dir.''')
        self.Button8.configure(command=lambda: openPostDir(self.Label4, self.Scrolledlistbox1))

        self.Button9 = Button(self.TNotebook1_pg2)
        self.Button9.place(relx=0.77, rely=0.54, height=94, width=87)
        self.Button9.configure(activebackground="#d9d9d9")
        self.Button9.configure(activeforeground="#000000")
        self.Button9.configure(background=_bgcolor)
        self.Button9.configure(disabledforeground="#a3a3a3")
        self.Button9.configure(foreground="#000000")
        self.Button9.configure(highlightbackground="#d9d9d9")
        self.Button9.configure(highlightcolor="black")
        self.Button9.configure(pady="0")
        self.Button9.configure(text='''Start plot''')
        self.Button9.configure(command=lambda: openPlot(self.Scrolledlistbox1))

        self.Label4 = Label(self.TNotebook1_pg2)
        self.Label4.place(relx=0.02, rely=0.04, height=21, width=310)
        self.Label4.configure(activebackground="#f9f9f9")
        self.Label4.configure(activeforeground="black")
        self.Label4.configure(background=_bgcolor)
        self.Label4.configure(disabledforeground="#a3a3a3")
        self.Label4.configure(foreground="#000000")
        self.Label4.configure(highlightbackground="#d9d9d9")
        self.Label4.configure(highlightcolor="black")
        self.Label4.configure(justify=LEFT)
        self.Label4.configure(relief=SUNKEN)
        self.Label4.configure(text='''Choose directory!''')

        self.Button12 = Button(self.TNotebook1_pg2)
        self.Button12.place(relx=0.77, rely=0.17, height=24, width=87)
        self.Button12.configure(activebackground="#d9d9d9")
        self.Button12.configure(activeforeground="#000000")
        self.Button12.configure(background=_bgcolor)
        self.Button12.configure(disabledforeground="#a3a3a3")
        self.Button12.configure(foreground="#000000")
        self.Button12.configure(highlightbackground="#d9d9d9")
        self.Button12.configure(highlightcolor="black")
        self.Button12.configure(pady="0")
        self.Button12.configure(text='''Reload dir.''')
        self.Button12.configure(command=lambda: reloadPostDir(self.Scrolledlistbox1))

        strDt = StringVar()
        self.Entry1 = Entry(self.TNotebook1_pg2)
        self.Entry1.place(relx=0.84, rely=0.3, relheight=0.09, relwidth=0.13)
        self.Entry1.configure(background="white")
        self.Entry1.configure(disabledforeground="#a3a3a3")
        self.Entry1.configure(font="TkFixedFont")
        self.Entry1.configure(foreground="#000000")
        self.Entry1.configure(insertbackground="black")
        self.Entry1.configure(width=54)
        self.Entry1.configure(textvariable=strDt)
        strDt.set("1e-5")

        self.Label6 = Label(self.TNotebook1_pg2)
        self.Label6.place(relx=0.77, rely=0.3, height=21, width=24)
        self.Label6.configure(background=_bgcolor)
        self.Label6.configure(disabledforeground="#a3a3a3")
        self.Label6.configure(foreground="#000000")
        self.Label6.configure(text='''dt =''')
        self.Label6.configure(width=24)

        self.Label7 = Label(self.TNotebook1_pg2)
        self.Label7.place(relx=0.77, rely=0.41, height=21, width=84)
        self.Label7.configure(background=_bgcolor)
        self.Label7.configure(disabledforeground="#a3a3a3")
        self.Label7.configure(foreground="#000000")
        self.Label7.configure(relief=GROOVE)
        self.Label7.configure(text='''t =''')
        self.Label7.configure(width=84)


        root.after(2000, self.update_clock)

    def update_clock(self):
        if solverProc is not None:
            if solverProc.poll() is not None:
                self.Label5.configure(text='''Simulation terminated.''')
            else:
                self.Label5.configure(text='''Simulation is running.''')
        root.after(2000, self.update_clock)


# The following code is added to facilitate the Scrolled widgets you specified.
class AutoScroll(object):
    '''Configure the scrollbars for a widget.'''

    def __init__(self, master):
        #  Rozen. Added the try-except clauses so that this class
        #  could be used for scrolled entry widget for which vertical
        #  scrolling is not supported. 5/7/14.
        try:
            vsb = ttk.Scrollbar(master, orient='vertical', command=self.yview)
        except:
            pass
        hsb = ttk.Scrollbar(master, orient='horizontal', command=self.xview)

        #self.configure(yscrollcommand=self._autoscroll(vsb),
        #    xscrollcommand=self._autoscroll(hsb))
        try:
            self.configure(yscrollcommand=self._autoscroll(vsb))
        except:
            pass
        self.configure(xscrollcommand=self._autoscroll(hsb))

        self.grid(column=0, row=0, sticky='nsew')
        try:
            vsb.grid(column=1, row=0, sticky='ns')
        except:
            pass
        hsb.grid(column=0, row=1, sticky='ew')

        master.grid_columnconfigure(0, weight=1)
        master.grid_rowconfigure(0, weight=1)

        # Copy geometry methods of master  (taken from ScrolledText.py)
        if py3:
            methods = Pack.__dict__.keys() | Grid.__dict__.keys() \
                  | Place.__dict__.keys()
        else:
            methods = Pack.__dict__.keys() + Grid.__dict__.keys() \
                  + Place.__dict__.keys()

        for meth in methods:
            if meth[0] != '_' and meth not in ('config', 'configure'):
                setattr(self, meth, getattr(master, meth))

    @staticmethod
    def _autoscroll(sbar):
        '''Hide and show scrollbar as needed.'''
        def wrapped(first, last):
            first, last = float(first), float(last)
            if first <= 0 and last >= 1:
                sbar.grid_remove()
            else:
                sbar.grid()
            sbar.set(first, last)
        return wrapped

    def __str__(self):
        return str(self.master)

def _create_container(func):
    '''Creates a ttk Frame with a given master, and use this new frame to
    place the scrollbars and the widget.'''
    def wrapped(cls, master, **kw):
        container = ttk.Frame(master)
        return func(cls, container, **kw)
    return wrapped

class ScrolledListBox(AutoScroll, Listbox):
    '''A standard Tkinter Text widget with scrollbars that will
    automatically show/hide as needed.'''
    @_create_container
    def __init__(self, master, **kw):
        Listbox.__init__(self, master, **kw)
        AutoScroll.__init__(self, master)

if __name__ == '__main__':
    vp_start_gui()

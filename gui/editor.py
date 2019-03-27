from Tkinter import *
from tkSimpleDialog import askstring
from tkFileDialog   import asksaveasfilename

from tkMessageBox import askokcancel

import tkFont

import os
import time


class Quitter(Frame):
    def __init__(self, parent=None, text=None):
        Frame.__init__(self, parent)
        self.pack()
        self.text = text
        widget = Button(self, text='Close editor', command=self.quit)
        widget.pack(expand=YES, fill=BOTH, side=LEFT)

    def quit(self):
        if (self.text.edit_modified()):
            ans = askokcancel('Really close the file?', "There were unsaved changes! Do you want to close without saving?")
            if ans: Frame.quit(self)
        else:
            Frame.quit(self)

class ScrolledText(Frame):
    fileChanged = False

    def __init__(self, parent=None, text='', file=None):
        Frame.__init__(self, parent)
        self.pack(expand=YES, fill=BOTH)
        self.makewidgets()
        self.settext(text, file)

    def makewidgets(self):
        sbar = Scrollbar(self)
        text = Text(self, relief=SUNKEN)
        sbar.config(command=text.yview)
        text.config(yscrollcommand=sbar.set)
        sbar.pack(side=RIGHT, fill=Y)
        text.pack(side=LEFT, expand=YES, fill=BOTH)
        self.text = text

    def settext(self, text='', file=None):
        if file:
            text = open(file, 'r').read()
        self.text.delete('1.0', END)
        self.text.insert('1.0', text)
        self.text.mark_set(INSERT, '1.0')
        self.text.edit_modified(False)
        self.text.focus()

    def gettext(self):
        return self.text.get('1.0', END+'-1c')


class SimpleEditor(ScrolledText):
    fileName=None
    statusText=None

    def __init__(self, parent=None, file=None):
        frm = Frame(parent)
        frm.pack(fill=X)
        Button(frm, text='Save',  command=self.onSave).pack(side=LEFT)
        Button(frm, text='Save as',  command=self.onSaveAs).pack(side=LEFT)
        Button(frm, text='Cut',   command=self.onCut).pack(side=LEFT)
        Button(frm, text='Paste', command=self.onPaste).pack(side=LEFT)
        Button(frm, text='Find',  command=self.onFind).pack(side=LEFT)        
        ScrolledText.__init__(self, parent, file=file)
        Quitter(frm, text=self.text).pack(side=LEFT)
        self.text.config(font=('FixedSys', 9, 'normal'))
        self.fileName = file
        self.statusText = StringVar()
        Label(frm, textvariable=self.statusText, bd=1, relief=SUNKEN, anchor=W).pack(side=RIGHT)
        self.statusText.set("Edit setup.")

    def onSave(self):
        if(self.text.edit_modified()):
            alltext = self.gettext()
            open(self.fileName, 'w').write(alltext)
            self.statusText.set("File saved at " + time.strftime("%H:%M:%S"))            
        self.text.edit_modified(False)

    def onSaveAs(self):
        options = {}
        options['defaultextension'] = '.ini'
        options['filetypes'] = [('ini files', '.ini')]
        options['initialdir'] = os.path.dirname(self.fileName)
        options['initialfile'] = os.path.basename(self.fileName)
        options['title'] = 'Save simulation setup'
        
        filename = asksaveasfilename(**options)

        if filename:
            alltext = self.gettext()
            open(filename, 'w').write(alltext)
            self.statusText.set("File saved at " + time.strftime("%H:%M:%S"))  
            self.fileName = filename
            self.text.edit_modified(False)

    def onCut(self):
        text = self.text.get(SEL_FIRST, SEL_LAST)
        self.text.delete(SEL_FIRST, SEL_LAST)
        self.clipboard_clear()
        self.clipboard_append(text)

    def onPaste(self):
        try:
            text = self.selection_get(selection='CLIPBOARD')
            self.text.insert(INSERT, text)
        except TclError:
            pass

    def onFind(self):
        target = askstring('SimpleEditor', 'Search String?')
        if target:
            where = self.text.search(target, INSERT, END)
            if where:
                print where
                pastit = where + ('+%dc' % len(target))
               #self.text.tag_remove(SEL, '1.0', END)
                self.text.tag_add(SEL, where, pastit)
                self.text.mark_set(INSERT, pastit)
                self.text.see(INSERT)
                self.text.focus()

if __name__ == '__main__':
    root = Tk()
    root.title("Configuration editor")
    root.geometry('{}x{}'.format(900, 400))
    root.wm_iconbitmap('@gui/main.xbm')
    try:
        SimpleEditor(file=sys.argv[1]).mainloop()
    except IndexError:
        SimpleEditor().mainloop()
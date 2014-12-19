'''
Created on 15.12.2014

@author: Lorenz
'''

if __name__ == '__main__':
    pass
#import Modules
import tkinter.filedialog
import os
import getpass
from tkinter import messagebox

#WINDOW-----------------------------------------------------------------------
#create a window
window = tkinter.Tk()
#name the window
window.title("AE Fit")
#set window standart size
window.geometry("600x600")
#set window icon
window.wm_iconbitmap('favicon.ico')

#VARIABLES----------------------------------------------------------------
user = getpass.getuser
var0 = 0
var1 = 0
var2 = 0
var3 = 0
var4 = 0

#ADD-FRAME-0--------------------------------------------------------------------
frame0 = tkinter.Frame()
frame0.pack(side = tkinter.LEFT)
frame0.place(x=5,y=5)

#ADD-FRAME-1---------------------------------------------------------------------
frame1 = tkinter.Frame()
frame1.pack(side = tkinter.BOTTOM)


#ADD-FRAME-2---------------------------------------------------------------------
frame2 = tkinter.Frame()
frame2.pack()
frame2.place(x=200,y=5)

#ADD-FRAME-3---------------------------------------------------------------------
frame3 = tkinter.Frame()
frame3.pack()
frame3.place(x=200,y=50)

#ADD-FRAME-4---------------------------------------------------------------------
frame4 = tkinter.Frame()
frame4.pack()
frame4.place(x=200,y=100)

#ADD-FRAME-5---------------------------------------------------------------------
frame5 = tkinter.Frame()
frame5.pack()
frame5.place(x=200,y=150)

#FILE-LIST-0---------------------------------------------------------------------
#define List
lis0 = tkinter.Listbox(frame0, yscrollcommand = "yes", height = 20, width = 30)
lis0.pack()

#FILE-LIST-1------------------------------------------------------------------------
lis1 = tkinter.Listbox(window, yscrollcommand = "yes", height = 20, width = 30)

#ADD-FILE-BUTTON-----------------------------------------------------------------
#define the browse function
def callback0():
    choosenpath = tkinter.filedialog.askopenfilename()
    choosenfile = os.path.split(choosenpath)[1]
    print(choosenpath)
    print(choosenfile)
    #writes file name in the list
    lis0.insert(1,choosenfile)
    lis1.insert(1,choosenpath)

#create a browse button
btn0 = tkinter.Button(frame0, text="add File", width = 11, command=callback0)
btn0.pack(padx=5,pady=5, side = tkinter.LEFT)

#DELET-FILE-BUTTON-------------------------------------------------------------
#define the delete function
def callback1():
    int1 = lis0.curselection()
    if not int1:
        messagebox.showwarning( "ERROR MESSAGE", "Sorry, but if you dosen't choose a file, I don't know what do do! ")
    else:
        lis0.delete(int1)
        lis1.delete(int1)
#create a delete file button
btn1 = tkinter.Button(frame0, text="delete File", width = 11, command=callback1)
btn1.pack(padx=5,pady=5, side = tkinter.LEFT)

#CLOSE-PROGRAMM-BUTTON--------------------------------------------------------------
#define the close function
def callback2():
    window.quit()
    
#create a close button
btn2 = tkinter.Button(frame1, text="close", width = 11, command=callback2)
btn2.pack(padx=5,pady=5, side = tkinter.LEFT)

#ENTRY-FIELD-0-------------------------------------------------------------------
def callback3(*event):
    var0 = ent0.get()
    lbl1.configure(text=var0)
    ent1.focus_set()
lbl0= tkinter.Label(frame2, text="Variable 1")
lbl1= tkinter.Label(frame2, text="")
ent0= tkinter.Entry(frame2)
btn3= tkinter.Button(frame2, text="set", command=callback3)

ent0.bind("<Return>", callback3)
    
lbl0.pack(padx=5,pady=5, side = tkinter.LEFT)
ent0.pack(padx=5,pady=5, side = tkinter.LEFT)
btn3.pack(padx=5,pady=5, side = tkinter.LEFT)
lbl1.pack(padx=5,pady=5, side = tkinter.LEFT)

#ENTRY-FIELD-1-------------------------------------------------------------------
def callback4(*event):
    var1 = ent1.get()
    lbl3.configure(text=var1)
    ent2.focus_set()
lbl2= tkinter.Label(frame3, text="Variable 2")
lbl3= tkinter.Label(frame3, text="")
ent1= tkinter.Entry(frame3)
btn4= tkinter.Button(frame3, text="set", command=callback4)

ent1.bind("<Return>", callback4)
    
lbl2.pack(padx=5,pady=5, side = tkinter.LEFT)
ent1.pack(padx=5,pady=5, side = tkinter.LEFT)
btn4.pack(padx=5,pady=5, side = tkinter.LEFT)
lbl3.pack(padx=5,pady=5, side = tkinter.LEFT)

#ENTRY-FIELD-1-------------------------------------------------------------------
def callback5(*event):
    var2 = ent2.get()
    lbl5.configure(text=var2)
    ent0.focus_set()
lbl4= tkinter.Label(frame4, text="Variable 3")
lbl5= tkinter.Label(frame4, text="   ")
ent2= tkinter.Entry(frame4)
btn5= tkinter.Button(frame4, text="set", command=callback5)

ent2.bind("<Return>", callback5)
    
lbl4.pack(padx=5,pady=5, side = tkinter.LEFT)
ent2.pack(padx=5,pady=5, side = tkinter.LEFT)
btn5.pack(padx=5,pady=5, side = tkinter.LEFT)
lbl5.pack(padx=5,pady=5, side = tkinter.LEFT)

#DRAW---------------------------------------------------------------------------
#draw the window
window.mainloop()
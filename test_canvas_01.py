from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
#from matplotlib.figure import Figure
#from Tkinter import *
import tkinter as Tk
#import numpy as np
#import math
#import matplotlib
#matplotlib.use('TkAgg')
#from matplotlib.backend_bases import key_press_handler
#import sys


root = Tk.Tk()

#Interface-----------------------------------------------------------

title_label = Tk.Button(root,text = "Add graph title", padx=2,pady=2)
xlabel = Tk.Button(root,text =      "Add X values    ", padx=2,pady=2)
ylabel = Tk.Button(root,text =      "Add Y values    ", padx=2,pady=2)
nameXaxis = Tk.Button(root,text =   "Name X axis      ", padx=2,pady=2)
nameYaxis = Tk.Button(root,text =   "Name Y axis      ", padx=2,pady=2)
meanLabel = Tk.Button(root,text = "Mean          ")
stderrorLabel = Tk.Button(root,text = "StdError:     ")

barGraph = Tk.Button(root,text = "Bar Graph  ",fg = "red", padx=2,pady=2)
lineGraph = Tk.Button(root,text = "Line Graph", fg = "red",padx=2,pady=2)
pieGraph = Tk.Button(root,text = "Pie Graph  ",fg = "red",padx=2,pady=2)

titleEntry = Tk.Entry(root)
xentry = Tk.Entry(root)
yentry = Tk.Entry(root)
nameXaxisEntry = Tk.Entry(root)
nameYaxisEntry = Tk.Entry(root)
meanText = Tk.Text(root,height=1,width=4)
stderrText = Tk.Text(root,height=1,width=4)        


title_label.grid(row = 0, column = 0,sticky = Tk.E)
xlabel.grid(row = 1, column = 0, sticky = Tk.E)
ylabel.grid(row = 2, column = 0, sticky = Tk.E)
nameXaxis.grid(row = 3, column = 0, sticky = Tk.E)
nameYaxis.grid(row = 4, column = 0, sticky = Tk.E)

barGraph.grid(row = 0,column = 1,ipadx=10,sticky=Tk.W)
lineGraph.grid(row = 1,column = 1,ipadx=10,sticky=Tk.W)
pieGraph.grid(row = 2,column = 1,ipadx=10,sticky=Tk.W)

meanLabel.grid(row = 3,column = 1,ipadx=10,sticky=Tk.W)
stderrorLabel.grid(row = 4,column = 1,ipadx=10,sticky=Tk.W)                      


titleEntry.grid(row = 0, column = 0,ipadx=100,sticky=Tk.W)
xentry.grid(row = 1, column = 0,ipadx=100,sticky=Tk.W)
yentry.grid(row = 2, column = 0,ipadx=100,sticky=Tk.W)
nameXaxisEntry.grid(row = 3, column = 0,ipadx=100,sticky=Tk.W)
nameYaxisEntry.grid(row = 4, column = 0,ipadx=100,sticky=Tk.W)
meanText.grid(row=3,column=3,sticky=Tk.W)
stderrText.grid(row=4,column=3,sticky=Tk.W)

# Adding line graph to Canvas--------------------------------------------


root.title("Naynts Graphs") 

#fig = Figure(figsize=(5,4), dpi=100)
fig = plt.figure(figsize=(5,4), dpi=100) 
ax = fig.add_subplot(111) 

canvas = FigureCanvasTkAgg(fig,root) 
#canvas.show() 
canvas.get_tk_widget().grid(row=7,column=0) 

toolbar_frame = Tk.Frame(root)
toolbar_frame.grid(row=9,column=0)
toolbar = NavigationToolbar2Tk(canvas, toolbar_frame) 
toolbar.update() 
#canvas._tkcanvas.grid(row=9,column=0)

# Adding features to graph

plt.xlabel('x label')
plt.ylabel('y label')
plt.title('Graph')
plt.show()

#root.mainloop()
Tk.mainloop()
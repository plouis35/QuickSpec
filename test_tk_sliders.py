# test_tk_sliders.py
import tkinter
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import numpy as np

# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure
from matplotlib.widgets import Slider, Button, RadioButtons, RangeSlider

root = tkinter.Tk()
root.wm_title("Embedding in Tk")
root.configure(background='grey')


def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

fig = Figure(figsize=(5, 4)) #, dpi=100)
hdul = fits.open("albireo-10.fit")
image_data = hdul[0].data
hdul.close()
fig, ax = plt.subplots()
im1 = ax.imshow(image_data, cmap='magma')
fig.set_facecolor("gray")
fig.colorbar(im1)
canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas.draw()

toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
#toolbar.config(background='grey')
#toolbar.config(foreground='gray')

#print(toolbar.toolitems)
toolbar.children['!button4'].pack_forget()
bt_quit = tkinter.Button(master=toolbar, text="Quit", command=_quit)
bt_quit.pack(side="left")

bt_load = tkinter.Button(master=toolbar, text="Load") #, command=_quit)
bt_load.pack(side="left")
toolbar.update()

min0 = 0
max0 = 25000
axmin = fig.add_axes([0.3, 0.9, 0.35, 0.03])
axmax = fig.add_axes([0.3, 0.95, 0.35, 0.03])

smin = Slider(axmin, 'Min : ', 0, 65535, valinit=min0, color = 'darkgray', orientation='horizontal')
smax = Slider(axmax, 'Max : ', 0, 65535, valinit=max0, color = 'darkgray', orientation='horizontal')

def update(val):
    im1.set_clim([smin.val,smax.val])
    fig.canvas.draw()

smin.on_changed(update)
smax.on_changed(update)

canvas.mpl_connect(
    "key_press_event", lambda event: print(f"you pressed {event.key}"))
canvas.mpl_connect("key_press_event", key_press_handler)

toolbar.pack(side=tkinter.TOP, fill=tkinter.X)
plt.subplots_adjust(bottom=0.2, top=0.9) 
canvas.get_tk_widget().pack(side=tkinter.BOTTOM, fill=tkinter.BOTH, expand=True)
tkinter.mainloop()


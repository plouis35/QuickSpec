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
from matplotlib.widgets import Slider, Button, RadioButtons

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
toolbar.update()

#print(toolbar.toolitems)
toolbar.children['!button4'].pack_forget()
bt_quit = tkinter.Button(master=toolbar, text="Quit", command=_quit)
bt_quit.pack(side="left")

bt_load = tkinter.Button(master=toolbar, text="Load") #, command=_quit)
bt_load.pack(side="left")

min0 = 0
max0 = 25000
axmin = fig.add_axes([0.2, 0.01, 0.35, 0.03])
axmax = fig.add_axes([0.2, 0.06, 0.35, 0.03])

#axmin = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
#axmax = fig.add_axes([0.15, 0.25, 0.0225, 0.63])


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
canvas.get_tk_widget().pack(side=tkinter.BOTTOM, fill=tkinter.BOTH, expand=True)
tkinter.mainloop()

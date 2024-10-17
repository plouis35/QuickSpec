import tkinter as tk  # PEP8: `import *` is not preferred
import tkinter.ttk as ttk
import tkinter.messagebox as msg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
from matplotlib.widgets import EllipseSelector, RectangleSelector


def select_callback(eclick, erelease):
    """
    Callback for line selection.

    *eclick* and *erelease* are the press and release events.
    """
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print(f"({x1:3.2f}, {y1:3.2f}) --> ({x2:3.2f}, {y2:3.2f})")
    print(f"The buttons you used were: {eclick.button} {erelease.button}")
    print(f"rotation={rect_selector.rotation}")


def toggle_selector(event):
    print('Key pressed.')
    if event.key == 't':
            if rect_selector.active:
                print(f'selector deactivated.')
                rect_selector.set_active(False)
            else:
                print(f'selector activated.')
                rect_selector.set_active(True)


class MyToolbar(NavigationToolbar2Tk):
    
    def __init__(self, canvas, parent):

        # list of toolitems to add/modify to the toolbar, format is:
        # (
        #   text, # the text of the button (often not visible to users)
        #   tooltip_text, # the tooltip shown on hover (where possible)
        #   image_file, # name of the image for the button (without the extension)
        #   name_of_method, # name of the method in NavigationToolbar2 to call
        # )
        self.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
            # change config button in toolbar
            ("Configure", "Change parameters", "subplots", "edit_config"),                    
            (None, None, None, None),
            ("Load", "Load a 2D spectrum ", "hand", "edit_config"),                    
            ('Save', 'Save the figure', 'filesave', 'save_figure'), 
        )
        super().__init__(canvas, parent) #, pack_toolbar=False)

    def edit_config(self):
        print("You have to create edit_parameters()")
        msg.showwarning("Warning", "You have to create edit_parameters()")
        
# --- main ---

y = [i**2 for i in range(101)]

root = tk.Tk()

fig = plt.Figure(figsize=(5, 5), dpi=100)
plot1 = fig.add_subplot(111)
plot1.plot(y)

props = dict(facecolor='blue', alpha=0.5)

rect_selector = RectangleSelector(
    plot1, select_callback,
    useblit=True,
    button=[1, 3], 
    minspanx=5, minspany=5,
    spancoords='pixels',
    interactive=True,
    props=props
    )
rect_selector.add_state('rotate')

fig.canvas.mpl_connect('key_press_event', toggle_selector)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
#canvas.get_tk_widget().grid(row=2, column=0)
canvas.get_tk_widget().grid(column=0, row=0, sticky=tk.NSEW)
root.grid_columnconfigure(0, weight=1)
root.grid_rowconfigure(0, weight=1)

#canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


###############    TOOLBAR    ###############
toolbarFrame = tk.Frame(master=root)
toolbarFrame.grid(row=0, column=0, sticky=tk.NW) #SEW)    
toolbar = MyToolbar(canvas, toolbarFrame)
#toolbar.grid(row=0, column=0, sticky=tk.N) #EW)


#toolbar.pack(side=tk.TOP, fill=tk.X)
#toolbar.update()
#canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

root.mainloop()

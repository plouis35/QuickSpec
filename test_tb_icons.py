from tkinter import *
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk, FigureCanvasTkAgg
from matplotlib.figure import Figure
#from matplotlib.backends._backend_tk import ToolTip

root = Tk()

class Navigator(NavigationToolbar2Tk):
    """
    Customized Navigator object
    - Removed mouse_move event.x and event.y display
    - Changed buttons layout
    """
    def mouse_move(self,event):
        pass
    
    """
    def _Button(self, text, file, command, extension='.gif'):
        im = PhotoImage(master=self, file=file)
        b = Button(master=self, text=text, padx=2,image=im, command=command,bg="DeepSkyBlue4",relief="flat")
        b._ntimage = im
        b.pack(side=LEFT)
        self.tool_buttons.append(b)
        return b
    """


    def __init__(self, canvas, root):
        self.tool_buttons = []
        self.toolbar_icons = ["icons/home.png", #provide paths of your own icons
                              "icons/backward.png",
                              "icons/forward.png",
                              None,
                              "icons/pan.png",
                              "icons/zoom.png",
                              "icons/config.png",
                              None,
                              "icons/save.png"]
        
        xmin, xmax = self.canvas.figure.bbox.intervalx
        height, width = 50, xmax-xmin
        Frame.__init__(self, master=self.window,
                          width=500, height=int(height))
        self.update()
        num = 0
        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self._Spacer()
            else:
                try:
                    button = self._Button(text=text, file=self.toolbar_icons[num],
                                          command=getattr(self, callback))
                    if tooltip_text is not None:
                        ToolTip.createToolTip(button, tooltip_text)
                except IndexError:
                    pass
            num+=1
        self.pack(side=TOP, fill=X)

    def destroy(self, *args):
        Frame.destroy(self)


f = Figure()

canvas = FigureCanvasTkAgg(f, root)
canvas.get_tk_widget().pack(fill=BOTH, expand=True)

toolbar = NavigationToolbar2Tk(canvas, root)

custom_toolbar  = Navigator(canvas, root)
custom_toolbar.config(background="DeepSkyBlue4")
root.mainloop()

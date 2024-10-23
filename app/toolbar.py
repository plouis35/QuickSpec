import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from img.image import spec2d

class CustomToolbar(NavigationToolbar2Tk):
    def edit_config(self):
        msg.showwarning("Warning", "You have to create edit_parameters()")

    def load_function(self):
        pass

    def param_function(self):
        pass

    def __init__(self, canvas, parent):
        # list of toolitems to add/modify to the toolbar, format is:
        # (
        #   text, # the text of the button (often not visible to users)
        #   tooltip_text, # the tooltip shown on hover (where possible)
        #   image_file, # name of the image for the button (without the extension)
        #   name_of_method, # name of the method in NavigationToolbar2 to call
        # )
        # this is enforced by MPL lib - sould use a DICT otherwize...
        self.toolitems = (
            ('Home', 'Reset zoom to original view', 'home', 'home'),
            ('Back', 'Back to previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
            # change config button in toolbar
            (None, None, None, None),
            #("Configure", "Change parameters", "subplots", "edit_config"),                    
#            ("Load", "Load a 2D spectrum ", "hand", "edit_config"),                    
            ('Save', 'Save the figure', 'filesave', 'save_figure'), 
        )

        super().__init__(canvas, parent)
        #self._message_label.config(background='gray20', foreground='white')
        self.add_custom_buttons()


    def add_custom_buttons(self):
        self.bt_load_img = tk.PhotoImage(file='./icons/load.png')        
        self.load_button = tk.Button(self, text="Load", command=spec2d.load_image, #) #,
                                     image=self.bt_load_img) #, bg='white', fg='white')
        self.load_button.pack(side=tk.LEFT) #, padx=2, pady=2)


        self.bt_conf_img = tk.PhotoImage(file='./icons/conf.png')
        self.param_button = tk.Button(self, text="Configuration", command=self.param_function, #) #, bg='gray20', fg='white')
                                     image=self.bt_conf_img) #, bg='gray20', fg='white')
        self.param_button.pack(side=tk.LEFT) #, padx=2, pady=2)

        self.bt_run_img = tk.PhotoImage(file='./icons/run.png')
        self.run_button = tk.Button(self, text="Run processing", command=self.param_function, #) #, bg='gray20', fg='white')
                                     image=self.bt_run_img) #, bg='gray20', fg='white')
        self.run_button.pack(side=tk.LEFT) #, padx=2, pady=2)

        self.bt_watch_img = tk.PhotoImage(file='./icons/eye.png')
        self.watch_button = tk.Button(self, text="Watch mode", command=self.param_function, #) #, bg='gray20', fg='white')
                                     image=self.bt_watch_img) #, bg='gray20', fg='white')
        self.watch_button.pack(side=tk.LEFT) #, padx=2, pady=2)

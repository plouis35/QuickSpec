import tkinter as tk
from tkinter import ttk

class MyCustomWidget(ttk.Frame):
    def __init__(self, master, title=None, *args, **kwargs):
        super().__init__( master, *args, **kwargs)

        self.title = title

        self.buid_widgets()
        return
    
    def buid_widgets(self):
        self.main_frm = ttk.Frame(self)
        self.main_frm.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        
        if self.title is not None:
            ttk.Label(self.main_frm, text=self.title).pack(side=tk.LEFT)

        self.ent_cod = ttk.Entry(self.main_frm, width=10)
        self.ent_cod.pack(side=tk.LEFT, fill=tk.Y)

        self.btn_show_dt = ttk.Button(self.main_frm, text='OK')
        self.btn_show_dt.pack(side=tk.LEFT)

        self.ent_txt = ttk.Entry(self.main_frm)
        self.ent_txt.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        return


if __name__== '__main__':
    app = tk.Tk()
    app.geometry('600x200')

    s = ttk.Style()
    s.configure('new.TFrame', background='#7AC5CD')

    frm_main = ttk.Frame(app, style='new.TFrame')
    frm_main.pack(fill=tk.BOTH, expand=1)
    frm_main.grid_columnconfigure(0, weight=1)
    frm_main.grid_rowconfigure(0, weight=1)

    frm_container = ttk.Frame(frm_main)
    frm_container.grid(row=0, column=0, sticky=tk.NSEW)
    frm_container.columnconfigure(0, weight=1)
    frm_container.columnconfigure(1, weight=100)

    ttk.Label(frm_container, text='Test0:').grid(row=0, column=0)
    ttk.Entry(frm_container).grid(row=0, column=1, sticky=tk.EW)

    wd_custom = MyCustomWidget(frm_container, title='Test1:')
    wd_custom.grid(row=1, column=0, columnspan=2, sticky=tk.EW)

    tk.Text(frm_container).grid(row=2, column=0, columnspan=2, sticky=tk.NSEW)
    app.mainloop()

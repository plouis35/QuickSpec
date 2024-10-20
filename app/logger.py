import logging
import tkinter as tk

class LogHandler(logging.Handler):
    def __init__(self, text_widget):
        logging.Handler.__init__(self)
        logging.Handler.setFormatter(self,
            logging.Formatter('%(asctime)s %(levelname)s - %(message)s',
            datefmt='%H:%M:%S'))
   
        self.text_widget = text_widget

    def emit(self, record):
        msg = self.format(record)
        
        def append():
            self.text_widget.configure(state='normal')
            self.text_widget.insert(tk.END, msg + '\n')
            self.text_widget.configure(state='disabled')
            self.text_widget.yview(tk.END)
            
        self.text_widget.after(0, append)

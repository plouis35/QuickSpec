import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import logging
import matplotlib.pyplot as plt

# Utiliser le style 'dark_background' de Matplotlib
plt.style.use('dark_background')

class CustomToolbar(NavigationToolbar2Tk):
    def __init__(self, canvas, parent):
        super().__init__(canvas, parent)
        self.config(background='gray20')
        self._message_label.config(background='gray20', foreground='white')
        self.add_custom_buttons()

    def add_custom_buttons(self):
        self.load_button = tk.Button(self, text="Load", command=self.load_function, bg='gray20', fg='white')
        self.load_button.pack(side=tk.LEFT, padx=2, pady=2)
        self.param_button = tk.Button(self, text="Paramètres", command=self.param_function, bg='gray20', fg='white')
        self.param_button.pack(side=tk.LEFT, padx=2, pady=2)

    def load_function(self):
        # Fonction vide pour le bouton 'Load'
        pass

    def param_function(self):
        # Fonction vide pour le bouton 'Paramètres'
        pass

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Application Tkinter avec Matplotlib")
        self.configure(bg='gray20')

        # Barre d'outils commune
        toolbar_frame = tk.Frame(self, bg='gray20')
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)

        # Figure Matplotlib avec deux sous-graphiques
        figure = Figure(figsize=(5, 4), dpi=100)
        plot1 = figure.add_subplot(211)
        data = np.random.rand(10, 10)
        plot1.imshow(data, cmap='viridis')
        plot1.set_title("Image avec imshow", color='white')

        plot2 = figure.add_subplot(212)
        plot2.plot([0, 1, 2, 3], [9, 4, 1, 0])
        plot2.set_title("Graphique linéaire", color='white')

        canvas = FigureCanvasTkAgg(figure, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = CustomToolbar(canvas, toolbar_frame)
        toolbar.update()

        # Zone de sortie des logs en bas
        log_frame = tk.Frame(self, bg='gray20')
        log_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        log_text = tk.Text(log_frame, state='disabled', height=10, bg='gray20', fg='white')
        log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar = ttk.Scrollbar(log_frame, command=log_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        log_text['yscrollcommand'] = scrollbar.set

        # Configuration du logger
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG)
        handler = TextHandler(log_text)
        self.logger.addHandler(handler)

        # Ajout de quelques messages de log pour tester
        self.logger.debug("Debug message")
        self.logger.info("Info message")
        self.logger.warning("Warning message")
        self.logger.error("Error message")
        self.logger.critical("Critical message")

class TextHandler(logging.Handler):
    def __init__(self, text_widget):
        logging.Handler.__init__(self)
        self.text_widget = text_widget

    def emit(self, record):
        msg = self.format(record)
        def append():
            self.text_widget.configure(state='normal')
            self.text_widget.insert(tk.END, msg + '\n')
            self.text_widget.configure(state='disabled')
            self.text_widget.yview(tk.END)
        self.text_widget.after(0, append)

if __name__ == "__main__":
    app = App()
    app.mainloop()

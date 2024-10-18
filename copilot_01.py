import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import logging

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Application Tkinter avec Matplotlib")

        # Barre d'icônes
        toolbar = tk.Frame(self, bd=1, relief=tk.RAISED)
        icon = tk.Button(toolbar, text="Icone 1")
        icon.pack(side=tk.LEFT, padx=2, pady=2)
        toolbar.pack(side=tk.TOP, fill=tk.X)

        # Image Matplotlib avec imshow
        figure1 = Figure(figsize=(5, 2), dpi=100)  # Taille réduite
        plot1 = figure1.add_subplot(111)
        data = np.random.rand(10, 10)
        plot1.imshow(data, cmap='viridis')
        canvas1 = FigureCanvasTkAgg(figure1, self)
        canvas1.draw()
        canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Graphique Matplotlib
        figure2 = Figure(figsize=(5, 2), dpi=100)  # Taille réduite
        plot2 = figure2.add_subplot(111)
        plot2.plot([0, 1, 2, 3], [9, 4, 1, 0])
        canvas2 = FigureCanvasTkAgg(figure2, self)
        canvas2.draw()
        canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Barre d'outils commune pour les deux figures
        toolbar_frame = tk.Frame(self)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        toolbar = NavigationToolbar2Tk(canvas1, toolbar_frame)
        toolbar.update()

        # Fonction pour changer la figure active
        def change_figure(event, canvas):
            toolbar.canvas = canvas
            toolbar.update()

        # Associer les événements de clic pour changer la figure active
        canvas1.mpl_connect("button_press_event", lambda event: change_figure(event, canvas1))
        canvas2.mpl_connect("button_press_event", lambda event: change_figure(event, canvas2))

        # Zone de sortie des logs en bas
        log_frame = tk.Frame(self)
        log_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        log_text = tk.Text(log_frame, state='disabled', height=10)
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

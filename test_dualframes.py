import tkinter as tk
from tkinter import ttk
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FITS Image and Spectrum Viewer")

        # Create frames
        self.top_frame = ttk.Frame(self)
        self.bottom_frame = ttk.Frame(self)
        self.top_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        # Display FITS image
        self.display_fits_image("albireo-10.fit")

        # Display Spectrum
        self.display_spectrum([1, 2, 3, 4], [10, 20, 30, 40])  # Example data

        # Display quit button
        self.button = ttk.Button(master=self, text="Quit", command=self._quit)
        self.button.pack(side=tk.TOP)


    def display_fits_image(self, fits_path):
        hdul = fits.open(fits_path)
        image_data = hdul[0].data
        hdul.close()

        self.fig, ax = plt.subplots()
        ax.imshow(image_data, cmap='gray')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.top_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Add toolbar for zoom functionality
        toolbar = NavigationToolbar2Tk(self.canvas, self.top_frame)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


    def display_spectrum(self, x_data, y_data):
        fig, ax = plt.subplots()
        ax.plot(x_data, y_data)
        canvas = FigureCanvasTkAgg(fig, master=self.bottom_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


    def _quit(self):
        self.quit()     # stops mainloop
        self.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

if __name__ == "__main__":
    app = App()
    app.mainloop()



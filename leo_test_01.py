import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import logging

# Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create a console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)

# Create a file handler
file_handler = logging.FileHandler('example.log')
file_handler.setLevel(logging.DEBUG)

# Create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Add the format to the handlers
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(console_handler)
logger.addHandler(file_handler)

# Create the main window
root = tk.Tk()

# Create a figure and a canvas
fig = Figure(figsize=(5, 4), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Create the toolbar
toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()

# Pack the toolbar at the top of the window
toolbar.pack(side=tk.TOP, fill=tk.X)

# Create a Text widget for the log
log_frame = tk.Frame(root, width=50, height=10)
log_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
log = tk.Text(log_frame)
log.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

# Add some log messages
logger.info('This is an info message')
logger.warning('This is a warning message')
logger.error('This is an error message')

# Run the main loop
root.mainloop()

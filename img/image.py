import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Axes

class spec2d:
    def __init__(self, ax: Axes):
        self._ax = ax
        dummy_img = np.random.rand(100, 100)
        self._ax.imshow(dummy_img, cmap='viridis')
    

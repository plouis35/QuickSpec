import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Axes

class spec1d:
    def __init__(self, ax: Axes):
        self._ax = ax
        self._ax.plot([0, 1, 2, 3], [9, 4, 1, 0])

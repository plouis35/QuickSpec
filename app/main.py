import logging
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.image import AxesImage
from matplotlib.widgets import Button

from app.config import Config
from img.img_utils import ImgTools
from spc.calibration import Calibration

class Application(object):
    def __init__(self) -> None:
        self.conf = Config()
        self.create_gui()

    def run(self) -> None:
        plt.show(block=True)

    def create_gui(self) -> None:
        plt.style.use(self.conf.get_str('window', 'theme'))

        # create a single figure for both image & spectrum horizontaly packed
        self.figure, (self.axe_img, self.axe_spc) = plt.subplots(2, 1, figsize=(8, 6), num="QuickSpec")
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.99)      
        logging.debug('figure created')

        # intialize imgTools
        ImgTools(self.axe_img, self.axe_spc)

        # create load button
        self._load_image_ax = self.figure.add_axes(rect=(0.08, 0.54, 0.06, 0.05)) # (left, bottom, width, height
        self._bt_load_image = Button(self._load_image_ax, 'Load', color='k') #, hovercolor='grey')
        self._bt_load_image.on_clicked(ImgTools.open_image)

        # prepare calibration
        _calib = Calibration(self.axe_img, self.axe_spc)
        
        # create run button
        self._run_calib_ax: Axes = self.figure.add_axes(rect=(0.15, 0.54, 0.06, 0.05)) # (left, bottom, width, height
        self._bt_run = Button(self._run_calib_ax, 'Run', color='k') #, hovercolor='grey')
        self._bt_run.on_clicked(_calib.do_calibration)

        logging.debug('buttons created')


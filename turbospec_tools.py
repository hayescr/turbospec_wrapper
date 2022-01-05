import os
import subprocess
import numpy as np
from spec_tools import Spectrum


class SpectrumTS(Spectrum):
    def __init__(self, filepath, normed=True):
        if normed:
            usecols = (0, 1)
        else:
            usecols = (0, 2)

        data = np.genfromtxt(filepath, dtype=None, encoding=None,
                             usecols=usecols)
        wave = data[:, 0]
        flux = data[:, 1]
        super().__init__(wave, flux)

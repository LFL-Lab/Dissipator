# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:35:31 2022

@author: Evangelos Vlachos <evlachos@usc.edu>

"""

import numpy as np
from datetime import date,datetime
from qm.qua import dual_demod
#%% res_demod
def res_demod(I, Q):

    return (dual_demod.full("integW_cos", "out1", "integW_minus_sin", "out2", I),
            dual_demod.full("integW_sin", "out1", "integW_cos", "out2", Q))

def IQ_imbalance(g, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1-g**2)*(2*c**2-1))
    return [float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]

def delayed_gauss(amp, length, sigma):
    gauss_arg = np.linspace(-sigma, sigma, length)
    delay = 16 - length - 4
    if delay < 0:
        return amp * np.exp(-(gauss_arg ** 2) / 2)

    return np.r_[np.zeros(delay), amp * np.exp(-(gauss_arg ** 2) / 2), np.zeros(4)]
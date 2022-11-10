# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 15:40:58 2021

@author: lfl
"""

import numpy as np
import matplotlib.pyplot as plt

def convert_V_to_dBm(data_V):
    
    return 10*np.log10(1e3*data_V**2/50)
    

    



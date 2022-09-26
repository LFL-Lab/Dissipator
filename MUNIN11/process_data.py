# -*- coding: utf-8 -*-

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Import data for fitting
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
from meas_utilities import *
from fitting_test import *

pathR1 = 'D:\\Users\\lfl\\data\\MUNIN11\\res_6.81686GHz\\'
pathR2 = 'D:\\Users\\lfl\\data\\MUNIN11\\res_7.02625GHz\\'
pathR3 = 'D:\\Users\\lfl\\data\\MUNIN11\\res_7.1069GHz\\'

datapath = pathR3

foldername = 'ramsey_08-17-22\\'
name = "n_avg_12000_t0_4_tf_4000.0_dt_8_resettime_50000.0_08-17-22_164457"

data = getdata(datapath + foldername + name + ".csv")

times = np.array(data["times"]) * 4 / 1e3; # times in microseconds
quad = "Q"
I = np.array(data[quad]); 

fit, fig = perform_fit(expCosine, times, I, plot=True, 
                  maketitle=True, freq = 1.0,
                  title = "T2 ramsey", precision = 4, 
                  xlabel = "t (us)", ylabel = quad  + " (V)")
# fit = perform_fit(times, I, fitfunc = expCosine, plot=True, freq_0 = 0.25);

fig.savefig(datapath + foldername + name + "_fit.png")

print(fit)

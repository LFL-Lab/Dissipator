# -*- coding: utf-8 -*-

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Import data for fitting
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt
from meas_utilities import *
from fitting import *

datapath = r"D:\\Users\\lfl\\data\\candle\\res_7.13135GHz\\"

foldername = 'echo_09-17-22\\'
name = "n_avg_2000_t0_0_tf_200000.0_dt_250_resettime_400000.0_09-17-22_021605"

data = getdata(datapath + foldername + name + ".csv")

times = np.array(data["times"]) * 4 / 1e3; # times in microseconds
quad = "I"
I = np.array(data[quad]); 

fit, fig = perform_fit(Exp, times, I, plot=True, 
                  maketitle=True,
                  title = "echo", precision = 4, 
                  xlabel = "t (us)", ylabel = quad  + " (V)")
# fit = perform_fit(times, I, fitfunc = expCosine, plot=True, freq_0 = 0.25);

fig.savefig(datapath + foldername + name + "_fit.png")

print(fit)

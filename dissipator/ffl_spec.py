# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:13:40 2023

@author: lfl
"""

from dissipator import *
import instrument_init as inst
import plot_functions as pf
import h5py
from datetime import datetime
import os
from instrument_init import init_sa, init_sa_by_serial_number
import matplotlib.pyplot as plt



qb = dissipator('diss09', device_name='diss09_5578')
qb.update_value('ffl_atten', 40)
# flux_list = [50e-6, 60e-6, 70e-6, 80e-6, 90e-6]
fmin = 2.4e9
fmax = 4.7e9
df = 10e6
freqs = np.arange(fmin, fmax + df/2, df)
Idata = np.zeros((len(flux_list), len(freqs)))
Qdata = np.zeros((len(flux_list), len(freqs)))
for i, flux in enumerate(tqdm(flux_list)):
    inst.set_flux_bias(flux)
    I, Q, freqs, job = qb.ffl_spec_sweep_lo(fmin = fmin,fmax = fmax,df = df) 
    Idata[i] = I
    Qdata[i] = Q
    
for i, flux in enumerate(flux_list):
    plt.plot(freqs, np.abs(Idata[i] + 1j* Qdata[i]), label=f'flux={round(flux*1e6)} uA')
plt.legend()
plt.show()
# if __name__ == "__main__":
#     main()
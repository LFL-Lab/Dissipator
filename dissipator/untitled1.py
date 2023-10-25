# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 10:13:23 2023

@author: lfl
"""
import h5py
import numpy as np
from dissipator import dissipator
from qubit import qubit
import plot_functions as pf
import matplotlib.pyplot as plt
qb = dissipator('diss08_11a',device_name='diss08_11a')
dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\'
filename = 'cavity_ringdown_flux_ffl_len=0.0mA_fflFreq=3.05GHz_DA=30dB_fDA=13dB_rrLen=1000clks_navg=100000_3.h5'
f = h5py.File(dataDir + filename, 'r')
f.keys()
#ffl_len_list=[500,1000,1800,2800,4000]
ffl_len_list=[4000]
T2vals=np.zeros((1,11))
errors=np.zeros((1,11))
for i, ffl_len in enumerate(list(f.keys())):
    
    for j, amp in enumerate(list(f[ffl_len].keys())):
        
        I=f[ffl_len][amp]['I'][()]
        Q=f[ffl_len][amp]['Q'][()]
        ydata = np.abs(I+1j*Q)
        t_arr=f[ffl_len][amp]['time'][()]
        fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='cavity-reset',dt=t_arr[-1]*1e-6/len(t_arr))
        
        #fig = pf.plot_data(t_arr,ydata,sequence='cavity-reset',fitted_pars=fitted_pars,nAverages=2000, pi2Width=qb.pars['pi_half_len'],
                         #qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=0,amp_ffl_scale=j*0.1, amp=0.5,flux=0.7, error=error, ffl_len=ffl_len_list[i])
        
        T2vals[i,j]=fitted_pars[1]
        errors[i,j]=error[1]
        

fig, ax = plt.subplots(figsize=(6,6))
plot=ax.imshow(T2vals, interpolation='nearest', extent=[0,1,500,4000], vmin=0.3, vmax=7, aspect="auto", origin="lower")
#ax.set_yticklabels([500,1000,1800,2800,4000])
plt.xlabel("ffl_amp")
plt.ylabel("ffl time (ns) before ramsey")
plt.title('flux=850uA, ffl_freq=2.84')
plt.colorbar(plot)

fig, ax = plt.subplots(figsize=(6,6))
plt.scatter(np.linspace(0,1,11), T2vals[0,:])
plt.errorbar(np.linspace(0,1,11), T2vals[0,:], yerr=errors[0,:], fmt="o")
plt.xlabel("ffl_amp")
plt.ylabel("T2 (us)")
plt.title('flux=1.15mA, ffl_freq=2.84, ffl len=4000')
ax.set_ylim([0.3,8])

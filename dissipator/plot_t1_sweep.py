# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 16:41:53 2023

@author: lfl
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import h5py
import plot_functions as pf

#hf=h5py.File(f"G:\\Shared drives\\CavityCooling\\data\\diss08_07A\\20230407\\T1wFFL\\T1wFFL_fflFreq=3d6918GHz_DA=35dB_navg=5000_1.h5",'r')
path=f'G:\\Shared drives\\CavityCooling\\data\\diss08_07A\\20230407\\T1wFFL\\'
name='T1wFFL_fflFreq=3d6918GHz_DA=35dB_navg=5000_1.h5'
hf=h5py.File(path + name,'r')
fig,axs=plt.subplots(4,3)
i=1
j=1
for dset in hf['sweep_ffl_amp_13:16:09'].keys():
    I=(hf['sweep_ffl_amp_13:16:09'][dset])['I'][:]
    Q=(hf['sweep_ffl_amp_13:16:09'][dset])['Q'][:]
    t=(hf['sweep_ffl_amp_13:16:09'][dset])['t'][:]
    ydata = np.abs(I+1j*Q)
    
    fitted_pars, error = pf.fit_data(t,ydata,sequence='T1',dt=t[-1]*1e-6/len(t))
    #fig = pf.plot_data(t,ydata,sequence='T1',fitted_pars=fitted_pars,nAverages=5000, pi2Width=qb.pars['pi_half_len'],
                       #qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
    axs[i,j].plot(t, ydata*1e3, '-o', markersize = 3, c='C0')
    #ax.set_ylabel('Digitizer Voltage (mV)')
    #ax.set_xlabel('Delay ($\mu$s)')
    axs[i,j].plot(t,pf.decay(t, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
    #textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.4f GHz\n$T_1$ = %.3f $\mu$s\n$\hat{n}$ = %d'%(pi2Width,qubitDriveFreq*1e-9,fitted_pars[1],nAverages)
    #axs[i].set_title('T1 Measurement %03d' %(iteration))
    
    
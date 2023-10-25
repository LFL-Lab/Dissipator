# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 14:08:54 2023

@author: lfl
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import scipy as scy
# from matplotlib import cm
import csv
import itertools
from scipy.interpolate import interp1d
import scipy.fftpack
import time
import os
from json import loads
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.colors import LightSource
from types import SimpleNamespace
pi=np.pi
import seaborn as sns; sns.set() # styling
from matplotlib.ticker import FormatStrFormatter
from scipy.signal import butter,lfilter,freqz,find_peaks,peak_widths
# import imageio

# set some deafault
# for a complete set of parameters "print(plt.rcParams)"
sns.set_style('ticks')
plt.rcParams['font.family'] =  'Arial'
plt.rcParams['font.size'] = 16
plt.rcParams['figure.dpi'] = 150
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.major.top"] = True
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams["xtick.major.bottom"] = True
plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.bottom"] = True
plt.rcParams["ytick.left"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams["ytick.major.right"] = True
plt.rcParams["ytick.labelright"] = False
plt.rcParams["ytick.minor.visible"] = False



def plot_data_inter(x_vector,y_b,y_c,y_cf,seq='interleaved_cavity_reset',qubitDriveFreq=3.8e9,fflDriveFreq=2e9,
                              nAverages=1,
                              stepSize=5e-6,fitted_pars_c=np.zeros(7),fitted_pars_cf=np.zeros(7),
                              fitted_pars_b=np.zeros(7), savefig=True, amp=1, ffl_atten=0, rr_atten=0, flux=0, amp_ffl_scale=0, error_b=[0,0,0,0,0],error_c=[0,0,0,0,0],error_cf=[0,0,0,0,0], ffl_len=0.):
    #x_vector = x_vector
    y_b = y_b*1e3
    y_c = y_c*1e3
    y_cf = y_cf*1e3
    #print(len(x_vector))
    #print(len(y_b))
    lim=min(y_b)
    
    #amp = (max(y_vector)-min(y_vector))/2
    #offset = np.mean(y_vector)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x_vector, y_b, marker='o', label='bare, T2=%.2f $\pm$ %.2f $\mu$s'%(fitted_pars_b[1], error_b[1]), color='r')
    ax.scatter(x_vector, y_c, marker='x', label='with cavity, T2=%.2f $\pm$ %.2f $\mu$s'%(fitted_pars_c[1], error_c[1]), color='b')
    ax.scatter(x_vector, y_cf, marker='+', label='with cavity and ffl, T2=%.2f $\pm$ %.2f $\mu$s'%(fitted_pars_cf[1], error_cf[1]), color='g')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Separation ($\mu$s)')
    ax.plot(x_vector,decay(x_vector, fitted_pars_b[0], fitted_pars_b[1], fitted_pars_b[2]),'r', color='r')
    ax.plot(x_vector,decay(x_vector, fitted_pars_c[0], fitted_pars_c[1], fitted_pars_c[2]),'r', color='b')
    ax.plot(x_vector,decay(x_vector, fitted_pars_cf[0], fitted_pars_cf[1], fitted_pars_cf[2]),'r', color='g')
    
    textstr = '$\omega_d$ = %.4f GHz\n$\hat{n}$ = %d \n$flux$= %.3f mA\n ffl scale=%.2f\n rr scale=%.2f \n $\omega_{ffl}$=%.3f GHz \n$(ffl,rr) atten$=(%.1f, %.1f) db \n ffl len=%.1f'%(qubitDriveFreq*1e-9,nAverages, flux,amp_ffl_scale, amp, fflDriveFreq*1e-9, ffl_atten,rr_atten, ffl_len)
    if seq=='interleaved_cavity_reset':
        ax.set_title('Interleaved Cavity Reset')
    else:
        ax.set_title('Interleaved Cavity Cooling')
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)
    plt.gca().set_ylim(bottom=lim)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return fig

def decay(x,amp,tau,offset):
    return amp*np.exp(-x/tau)+offset
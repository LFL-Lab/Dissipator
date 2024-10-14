# -*- coding: utf-8 -*-
"""
Created on Mon May 24 13:13:39 2021

@author: Evangelos Vlachos <evlachos@usc.edu>
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
from Utilities import convert_V_to_dBm
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

#%% power_plot
def power_plot(freqs,signal,power,fc):
    plt.plot(freqs*1e-6, signal,'-')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Power [dBm]')
    plt.show()

#%% tof_plot
def tof_plot(adc1,adc2,delay=0,offsets=[0,0]):
    plt.figure()
    plt.title('time-of-flight calibration analysis')
    plt.plot(adc1)
    plt.gcf().text(1, 0.15, f"Electrical delay: {delay}\nChannel 1 offset: {offsets[0]*1e3:.1f} mV\nChannel 2 offset: {offsets[1]*1e3:.1f} mV", fontsize=14)
    plt.plot(adc2)
    plt.legend(["adc1", "adc2"])
    plt.show()

#%% spec_plot
def resonator_spec_plot(data,qb_pars,fwhm=0,fc=0,iteration=1,**kwargs):
    freq = data['freqs']*1e-9
    df = (freq[1]-freq[0])*1e9
    I = data['I']*1e3
    Q = data['Q']*1e3
    mag = np.abs(I+1j*Q)
    power=10*np.log10((10**-3)*(mag**2)/50)

    phase = np.unwrap(np.angle(I+1j*Q,deg=True),period=360)
    fig = plt.figure(figsize=(8,5))

# Power data
    ax1 = fig.add_subplot(221)
    ax1.plot(freq, power, '-o', markersize=3, c='C0')
    ax1.set_xlabel('Frequency (GHz)')
    ax1.set_ylabel('Magnitude (dBm)')

    # Phase data
    ax2 = fig.add_subplot(222)
    ax2.plot(freq, phase, '-o', markersize=3, c='C0')
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Phase (deg)')

    # Additional subplot on the bottom row
    
    ax3 = fig.add_subplot(223)
    ax3.plot(freq, I, '-o', markersize=3, c='r', label='I')
    ax3.set_xlabel('Frequency (GHz)')
    ax3.set_ylabel('Voltage (mV)')
    ax3.set_title('Spectroscopy Data (I)')
    ax4 = fig.add_subplot(224)
    ax4.plot(freq,Q , '-o', markersize=3, c='b', label='Q')
    ax4.set_xlabel('Frequency (GHz)')
    ax4.set_ylabel('Voltage (mV)')
    ax4.set_title('Spectroscopy Data (Q)')
    

    txt = f'$\omega_c$ = {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-3:.3f} kHz\n$\kappa$ = {2*np.pi*fwhm*1e-6:.3f} MHz\nReadout attenuation: {qb_pars["readout_atten"]} dB\ndf = {df*1e-3:.1f} kHz'
    
    plt.gcf().text(1, 0.15, txt, fontsize=14)
    # fig.set_title(f'{element} spectroscopy {iteration}')
    plt.tight_layout()

def qubit_spec_plot(data,qb_pars,iteration=1,find_peaks=True, **kwargs):

    amp_q_scaling = kwargs['amp_q_scaling']
    rr_atten = kwargs['rr_atten']
    qubit_LO = kwargs['qubit_LO']
    width = kwargs['peak_width']

    freq = data['freqs']*1e-9
    df = freq[1]-freq[0]
    I = data['I']*1e3
    Q = data['Q']*1e3
    mag = np.abs(I+1j*Q)

    sigma = np.std(mag)
    print(f'Peak threshold at {np.mean(mag)+2*sigma}')
    peaks,_ = scy.signal.find_peaks(mag,height=np.mean(mag)+2*sigma,distance=200,width=width)
    try:
        for i in peaks:
            print(f'Peaks at: {round(freq[i],5)} GHz\n')
    except:
        print('Peaks not found or do not exist.')

    fig = plt.figure(figsize=(8,8))

# Power data
    ax1 = fig.add_subplot(2,2,(1,2))
    ax1.plot(freq, mag, '-o', markersize=3, c='C0')
    ax1.set_xlabel('Frequency (GHz)')
    ax1.set_ylabel('Magnitude (mV)')

    # Additional subplot on the bottom row
    ax2 = fig.add_subplot(223)
    ax2.plot(freq, I, '-o', markersize=3, c='r', label='I')
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Voltage (mV)')

    ax3 = fig.add_subplot(224)
    ax3.plot(freq, Q, '-o', markersize=3, c='r', label='Q')
    ax3.set_xlabel('Frequency (GHz)')
    ax3.set_ylabel('Voltage (mV)')

    if len(peaks) == 2:
        alpha = 2*(freq[peaks[0]]-freq[peaks[1]])*1e3
        txt = '$f_{01}$ = %.4f GHz\n$f_{02}$/2 = %.4f GHz\n$\\alpha$ = %.1f MHz\n$P_r$ = %.1f dBm\n$f_r$ = %.4f GHz\n$f_{LO}$=%.5f GHz\n$a_q$ = %.5f'%(freq[peaks[1]],freq[peaks[0]],alpha,rr_atten,qb_pars['rr_freq']*1e-9, qubit_LO*1e-9, amp_q_scaling)
    elif len(peaks) == 1:
        txt = '$f_{01}$ = %.4f GHz\n$P_r$ = %.1f dBm\n$f_r$ = %.4f GHz\n$f_{LO}$=%.5f GHz\n$a_q$ = %.5f'%(freq[peaks[0]], rr_atten, qb_pars['rr_freq']*1e-9, qubit_LO*1e-9,amp_q_scaling)
    else:
        txt = '$P_r$ = %.1f dBm\n$f_r$ = %.4f GHz\n$a_q$ = %.5f\n$f_{LO}$=%.5f GHz'%(rr_atten,qb_pars['rr_freq']*1e-9, amp_q_scaling, qubit_LO*1e-9)
    
    plt.gcf().text(1, 0.15, txt, fontsize=14)
    # fig.set_title(f'{element} spectroscopy {iteration}')
    plt.tight_layout()
    plt.show()

#%% init_IQ_plot
def init_IQ_plot():
    '''initialize axes for continuous plotting on the IQ plane'''
    plot = sns.jointplot()
    plot.set_axis_labels('I [mV]', 'Q [mV]')
    plot.ax_marg_x.grid('off')
    plot.ax_marg_y.grid('off')
    plot.fig.tight_layout()
    ax = plt.gca()
    return plot, ax

#%% heatmap_plot
def heatplot(xdata, ydata, data, xlabel = "", ylabel = "", normalize=False, cbar_label = 'log mag', **kwargs):
    fig,ax = plt.subplots(1,1,figsize=(4,3), dpi=300)
    if normalize:
        cbar_label += ' (normalized)'

    df = pd.DataFrame(data, columns = xdata, index = ydata)

    if normalize:
        df = df.apply(lambda x: (x-x.mean())/x.std(), axis = 1)

    cbar_options = {
        'label':    cbar_label,
        'ticks':    np.around(np.linspace(np.amin(data),np.amax(data),5),0),
        'pad':      0.1,
        'values':   np.linspace(np.amin(data),np.amax(data),1000),
        'shrink':   1.2,
        'location': 'right',
    }

    heatmap_opts = {
        'ax':    ax,
        'linewidths':  0,
        # 'xticklabels': np.linspace(min(xdata),max(xdata)+0.5,5),
        # 'yticklabels': np.linspace(min(ydata),max(ydata)+0.5,5),
        'vmin':        np.amin(data),
        'vmax':        np.amax(data)
    }

    hm = sns.heatmap(df, ax=ax,cmap = 'viridis', cbar_kws=cbar_options, **kwargs)
    hm.set_xlabel(xlabel, fontsize=12)
    hm.set_ylabel(ylabel, fontsize=12)
    hm.spines[:].set_visible(True)
    ax.tick_params(direction='out',length=0.01,width=0.5,bottom=True, top=True, left=True, right=True,labeltop=False, labelbottom=True,labelrotation=90,labelsize=8,size=8)
    plt.yticks(rotation=0)
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

def punchout_plot(data, xlabel = "", ylabel = "", 
             cbar_label = 'log mag', **kwargs):
    
    freqs = data['freqs']
    I = data['I']
    Q = data['Q']
    z_data = convert_V_to_dBm(np.array(data['mag']))
    phase = np.unwrap(np.angle(data['z_data'],deg=True),period=360)
    attenuations = np.array(data['attenuations'], dtype=int)
    chi= freqs[0][np.argmin(z_data[0])] - freqs[0][np.argmin(z_data[-1])]

    print(f'Dispersive shift: {0.5*chi/np.pi*1e9:.1f} kHz')

    fig = plt.figure(figsize=(6,5), dpi=300)
    ax1  = fig.add_subplot(6,1,(1,4))
    df = pd.DataFrame(z_data, columns = freqs[0], index = attenuations)
    
    cbar_options = {
        'label':    cbar_label,
        # 'ticks':    np.around(np.linspace(np.amin(z_data),np.amax(z_data),5),1),
        'pad':      0.05,
        # 'values':   np.linspace(np.amin(z_data),np.amax(z_data),1000),
        'shrink':   0.6,
        'location': 'top',

    }
  
    hm = sns.heatmap(df, ax=ax1,cmap = 'seismic', cbar_kws=cbar_options)
    hm.set_ylabel(ylabel, fontsize=10)
    hm.spines[:].set_visible(True)
    ax1.tick_params(direction='out',length=0.01,width=0.5,bottom=False, top=False, left=True, right=False,labeltop=False, labelbottom=True,labelrotation=90,labelsize=10,size=8)
    plt.yticks(rotation=0)
    
    ax2 = fig.add_subplot(6,1,(5,6))
    ax2.plot(freqs[0],z_data[0],'-', markersize = 2, c='b',label=f'$P_r$ = -{attenuations[0]} dB')
    ax3 = ax2.twinx()
    ax3.plot(freqs[0],z_data[-1],'-', markersize = 2, c='r',label=f'$P_r$ = -{attenuations[-1]} dB')
    ax2.legend()
    ax3.legend()
    ax2.set_xlabel(xlabel, fontsize=10)
    ax2.set_ylabel(cbar_label,fontsize=10)
    ax2.tick_params(axis='x',labeltop=False, labelbottom=True,labelrotation=90,labelsize=8,size=8)
    ax2.tick_params(axis='y',labeltop=False, labelbottom=True,labelrotation=0,labelsize=6,size=8)
    ax3.tick_params(axis='y',labeltop=False, labelbottom=True,labelrotation=0,labelsize=6,size=8)
    
    fc1,_ = fit_res(freqs[0],  I=I[0], Q=Q[0])
    fc2,_ = fit_res(freqs[0], I=I[-1], Q = Q[-1])
    txt = f'$f_1$ = {fc1:.6f} GHz\n$f_2$ = {fc2:.6f} GHz\n$2\chi/2\pi$ = {(fc2-fc1)*1e3:.3f} MHz'
    plt.gcf().text(0.95, 0.15, txt, fontsize=12)
    
    # plt.tight_layout()
        
    # return df


    # return hm;

#%% plot_single_shot
def plot_single_shot(datadict):
    plot, ax = init_IQ_plot()
    datadict = {key: np.array(value, dtype=float) for key,value in datadict.items()}
    datadict = {key: value*1e3 for key,value in datadict.items()} # convert to mV
    datadict = {key: value.tolist() for key,value in datadict.items()} # convert to list

    states = []
    # for key,value in datadict.items():
    #     print(key+':'+str(len(value))+'\n')
    [states.append(r'$|g\rangle$') for i in range(len(datadict['I']))]
    [states.append(r'$|e\rangle$') for i in range(len(datadict['Iexc']))]
    data = {
            'I [mV]':   np.hstack((datadict['I'],datadict['Iexc'])),
            'Q [mV]':   np.hstack((datadict['Q'],datadict['Qexc'])),
            'States':   states
                }
    I = np.array(datadict["I"])
    Q = np.array(datadict["Q"])
    Iexc=np.array(datadict["Iexc"])
    Qexc = np.array(datadict["Qexc"])
    
    y_gr = np.average(np.abs(I+1j*Q))
    phase_gr=np.average(np.arctan(Q/I))
    y_exc= np.average(np.abs(Iexc+1j*Qexc))
    phase_exc=np.average(np.arctan(Qexc/Iexc))
    # print("ground voltage:", y_gr)
    # print("excited volgate:", y_exc )
    # print("ground phase:", phase_gr)
    # print("excited phase:", phase_exc )
    print('contrast_mag',y_gr-y_exc)
    print('rotation',(phase_gr+phase_exc)/2)
    print("contrast_Q",y_gr*(np.abs(np.sin((phase_gr-phase_exc)/2)))+(y_exc*np.abs((np.sin((phase_gr-phase_exc))/2))))
    #print("contrast_I",y_gr*(np.abs(np.cos((phase_gr-phase_exc)/2)))-(y_exc*np.abs((np.cos((phase_gr-phase_exc))/2))))
    
          
    dataF = pd.DataFrame(data=data)
    plot = sns.jointplot(data=dataF, x='I [mV]',y='Q [mV]',hue='States',ax=ax,space=0)
    plt.show()

#%% plot_mixer_opt
def plot_mixer_opt(par1,par2,power_data,cal='LO',element='qubit',fc=5e9):
    par1 = np.around(par1*1e3,1)
    par2 = np.around(par2*1e3,1)

    par1 = par1.tolist()
    par2 = par2.tolist()
    df = pd.DataFrame(data=power_data,index=par1,columns=par2)

    hm = sns.heatmap(df,cbar_kws={'label': "Power [dBm]"})

    if cal == 'LO':
        hm.set_ylabel('I [mV]')
        hm.set_xlabel('Q [mV]')
    elif cal == 'SB':
        hm.set_ylabel('Gain Imbalance[x 1e-3]')
        hm.set_xlabel('Phase Imbalance[x 1e-3]')

    hm.spines[:].set_visible(True)
    hm.tick_params(direction='out',length=0.01,width=0.5,bottom=True, top=False, left=True, right=True,labeltop=False, labelbottom=True,labelrotation=90,labelsize=10,size=10)
    plt.yticks(rotation=0)
    plt.tight_layout()
    if element == 'qubit':
        plt.title(f'Qubit Mixer {cal} Calibration at {round(fc*1e-9,4)} GHz')
    elif element == 'rr':
        plt.title(f'Readout Mixer {cal} Calibration at {round(fc*1e-9,4)} GHz')
    plt.show()


#%% fit_res
def fit_res(f_data,I,Q,res_type='notch'):
    
    z_data = np.abs(I+1j*Q)
    fc = f_data[np.argmin(z_data)]
    if res_type == 'notch':
        z_data = -z_data-min(-z_data)
        idx = np.argwhere(np.diff(np.sign(z_data - 0.5*max(z_data)))).flatten()
        if len(idx) >=2:
            fwhm = f_data[idx[1]] - f_data[idx[0]]
        else:
            fwhm = np.nan
    print(f'Resonant Frequency: {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-3} kHz\nkappa = {2*np.pi*fwhm*1e-3:.3f} kHz')
    return fc,fwhm

# def find_peak(f_data,z_data):
    
#     return fc

#%% fit_data
def fit_data(x_vector, y_vector, sequence='rabi', verbose=0, **kwargs):
    '''
    fit experimental data

    sequence:          'Rabi','ramsey', 'T1' or 'T2'
    x_vector:           time data
    y_vector:           voltage data
    dt:                 sequence stepsize. Used for extracting the frequency of the data
    '''
    x_vector = x_vector
    y_vector = y_vector * 1e3

    amp = (max(y_vector) - min(y_vector)) / 2
    offset = np.mean(y_vector)

    if sequence == "t-rabi" or sequence == "qubit_temp" or sequence == 'qb-reset':
        return fit_rabi(x_vector, y_vector, amp, offset, kwargs['dt'], verbose)
    elif sequence == "p-rabi":
        return fit_p_rabi(x_vector, y_vector, amp, offset, kwargs['dt'], verbose)
    elif sequence == "ramsey" or sequence == "cavity-cooling-ramsey" or sequence == "ramsey_chi":
        return fit_ramsey(x_vector, y_vector, amp, offset, kwargs['dt'], kwargs['fitFunc'], verbose)
    elif sequence == "echo" or sequence == 'cavity-reset' or sequence == 'cavity-cooling':
        return fit_echo(x_vector, y_vector, amp, offset, verbose)
    elif sequence == "T1" or sequence == "T1diss":
        return fit_T1(x_vector, y_vector, amp, offset, verbose)
    elif sequence == 'ringdown_off' or sequence == 'ringdown_on':
        return fit_ringdown(x_vector, y_vector, amp, offset, verbose)
    else:
        raise ValueError("Unknown sequence type")

def fit_rabi(x_vector, y_vector, amp, offset, dt, verbose):
    fitFunction = rabi
    period = 1e3 / (extract_freq(x_vector * 1e3, y_vector, dt, plot=0))
    print('Period Initial Guess: %.1f ns' % (period))
    phase = pi
    x_vector = x_vector * 1e3
    lb = [0.1 * amp, 0.1 * period, 0, -2 * abs(offset)]
    ub = [10 * amp, 10 * period, 2 * pi, 2 * abs(offset)]
    p0 = [amp, period, phase, offset]
    return fit_curve(fitFunction, x_vector, y_vector, p0, lb, ub, verbose)

def fit_p_rabi(x_vector, y_vector, amp, offset, dt, verbose):
    fitFunction = rabi
    period = 1 / (extract_freq(x_vector, y_vector, dt , plot=0))
    print('Amplitude Initial Guess: %.3f' % (period))
    phase = pi
    lb = [0.1 * amp, 0.1 * period, 0, -2 * abs(offset)]
    ub = [10 * amp, 10 * period, 2 * pi, 2 * abs(offset)]
    p0 = [amp, period, phase, offset]
    return fit_curve(fitFunction, x_vector, y_vector, p0, lb, ub, verbose)

def fit_ramsey(x_vector, y_vector, amp, offset, dt, fitFunc, verbose):
    f = extract_freq(x_vector, y_vector, dt, plot=0)
    if x_vector[-1] > 20 and x_vector[-1] < 80:
        tau = 30
    elif x_vector[-1] > 80:
        tau = 80
    else:
        tau = 20
    phi = 0
    amp = abs(amp)
    if fitFunc != 'envelope':
        p0 = [amp, f, phi, tau, offset]
        lb = [0.75 * amp, 0, -pi, 0.01, -2 * abs(offset)]
        ub = [2 * amp, 2 * f, pi, 100, 2 * abs(offset)]
        fitFunction = ramsey
    else:
        tau = 1
        env = get_envelope(y_vector, dt, distance=100)
        env = env(x_vector) + offset
        p0 = [amp, tau, offset]
        if offset < 0:
            lb = [0.95 * amp, 0.1, 2 * offset]
            ub = [2 * amp, 15, 0.5 * offset]
        else:
            lb = [0.9 * amp, 0.1, 0.9 * offset]
            ub = [1.1 * amp, 15, 1.1 * offset]
        fitFunction = decay
        y_vector = env
    return fit_curve(fitFunction, x_vector, y_vector, p0, lb, ub, verbose)

def fit_echo(x_vector, y_vector, amp, offset, verbose):
    if x_vector[-1] < 10:
        tau = 2
        tau_ub = 20
    else:
        tau = 20
        tau_ub = 300
    amp = y_vector[0] - y_vector[-1]
    p0 = [amp, tau, offset]
    amp_bounds = [0.95 * amp, 1.05 * amp]
    off_bounds = [0.95 * offset, 1.05 * offset]
    lb = [min(amp_bounds), 0.1, min(off_bounds)]
    ub = [max(amp_bounds), tau_ub, max(off_bounds)]
    fitFunction = decay
    return fit_curve(fitFunction, x_vector, y_vector, p0, lb, ub, verbose)

def fit_T1(x_vector, y_vector, amp, offset, verbose):
    tau = 2
    amp = y_vector[0] - y_vector[-1]
    offset = y_vector[-1]
    if amp < 0:
        p0 = [amp, tau, offset]
        lb = [10 * amp, 0.1, -2 * abs(offset)]
        ub = [0.5 * amp, 300, 2 * abs(offset)]
    else:
        p0 = [amp, tau, offset]
        lb = [0.5 * amp, 0.1, -2 * abs(offset)]
        ub = [10 * amp, 300, 2 * abs(offset)]
    fitFunction = decay
    return fit_curve(fitFunction, x_vector, y_vector, p0, lb, ub, verbose)

def fit_ringdown(x_vector, y_vector, amp, offset, verbose):
    tau = 0.2
    amp = y_vector[0] - y_vector[-1]
    offset = y_vector[-1]
    p0 = [amp, tau, offset]
    fitFunction = decay
    return fit_curve(fitFunction, x_vector, y_vector, p0, verbose=verbose)

def fit_curve(fitFunction, x_vector, y_vector, p0, lb=None, ub=None, verbose=0):
    if lb is not None and ub is not None:
        try:
            fitted_pars, covar = sp.optimize.curve_fit(fitFunction, x_vector, y_vector, p0=p0, method='dogbox', bounds=[lb, ub], xtol=1e-12, maxfev=4e3)
            error = np.sqrt(abs(np.diag(covar)))
        except:
            print('Fit failed. ')
            fitted_pars = p0
            error = np.zeros(len(p0))
    else:
        try:
            fitted_pars, covar = sp.optimize.curve_fit(fitFunction, x_vector, y_vector, p0=p0, method='dogbox', xtol=1e-12, maxfev=40e3)
            error = np.sqrt(abs(np.diag(covar)))
        except:
            print('Fit failed. ')
            fitted_pars = p0
            error = np.zeros(len(p0))

    

    if verbose == 1:
        print('-' * 100)
        print('Initial Guess:', np.around(p0, 1))
        if lb is not None and ub is not None:
            print('Lower Bounds:', np.around(lb, 1))
            print('Upper Bounds:', np.around(ub, 1))
        print('Best Fit Pars:', np.around(fitted_pars, 1))
        print('Error:', np.around(error, 1))
        print('-' * 100)
    return fitted_pars, error


#%% plot_data
def plot_data(x_vector, y_vector, sequence='rabi', **kwargs):
    # y_vector = y_vector * 1e3

    if sequence == "p-rabi":
        return plot_p_rabi(x_vector, y_vector, **kwargs)
    elif sequence == "t-rabi" or sequence == "qubit_temp" or sequence == 'qb-reset':
        return plot_rabi(x_vector, y_vector, **kwargs)
    elif sequence == "ramsey" or sequence == "cavity-cooling-ramsey" or sequence == "ramsey_chi":
        return plot_ramsey(x_vector, y_vector, **kwargs)
    elif sequence == "echo":
        return plot_echo(x_vector, y_vector, **kwargs)
    elif sequence == 'cavity-reset':
        return plot_cavity_reset(x_vector, y_vector, **kwargs)
    elif sequence == 'cavity-cooling':
        return plot_cavity_cooling(x_vector, y_vector, **kwargs)
    elif sequence == 'ringdown_off' or sequence == 'ringdown_on':
        return plot_ringdown(x_vector, y_vector, **kwargs)
    elif sequence == "T1" or sequence == "dissT1":
        return plot_T1(x_vector, y_vector, **kwargs)
    else:
        raise ValueError("Unknown sequence type")
    
def plot_p_rabi(x_vector, y_vector, **kwargs):
    fig = plt.figure(figsize=(5,5))
    # I data
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(x_vector, y_vector[0]*1e3, '-o', markersize=3, c='C0')
    ax1.set_ylabel('I (mV)', fontsize=10)
    ax1.set_xlabel('Pulse Amplitude Scaling', fontsize=10)
    ax1.set_title('I Data', fontsize=10)
    ax1.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=4, labelsize=8)

    # Q data
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(x_vector, y_vector[1]*1e3, '-o', markersize=3, c='C1')
    ax2.set_ylabel('Q (mV)', fontsize=10)
    ax2.set_xlabel('Pulse Amplitude Scaling', fontsize=10)
    ax2.set_title('Q Data', fontsize=10)
    ax2.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=4, labelsize=8)

    # Magnitude data
    magnitude = np.sqrt(np.array(y_vector[0])**2 + np.array(y_vector[1])**2)
    ax3 = fig.add_subplot(2, 2, (3,4))
    ax3.plot(x_vector, magnitude*1e3, '-o', markersize=3, c='C2')
    ax3.set_ylabel('Magnitude (mV)', fontsize=10)
    ax3.set_xlabel('Pulse Amplitude Scaling', fontsize=10)
    ax3.set_title('Magnitude Data', fontsize=10)
    ax3.plot(x_vector, rabi(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2], kwargs['fitted_pars'][3]), 'r')
    ax3.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=4, labelsize=8)
    textstr = format_textstr_p_rabi(kwargs)
    plt.gcf().text(0.96, 0.45, textstr, fontsize=10)
    plt.tight_layout()
    # plt.show()


def plot_rabi(x_vector, y_vector, **kwargs):
    fig, ax = plt.subplots()
    ax.plot(x_vector * 1e3, y_vector*1e3, '-o', markersize=3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Duration (ns)')
    ax.plot(x_vector * 1e3, rabi(x_vector * 1e3, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2], kwargs['fitted_pars'][3]), 'r')
    ax.set_title('Rabi Measurement %03d' % (kwargs['iteration']))
    textstr = format_textstr_rabi(kwargs)
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)
    plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=8)


def plot_ramsey(x_vector, y_vector, **kwargs):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    fontSize = 16
    tickSize = 11
    markersize = 4
    linewidth = 2
    ax1.plot(x_vector, y_vector, 'o', markersize=markersize, c='C0')
    ax1.set_ylabel('Digitizer Voltage (mV)', fontsize=fontSize)
    ax1.set_xlabel('Pulse Separation ($\mu$s)')
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        label.set_fontsize(tickSize)
    if kwargs['fitFunc'] == 'envelope':
        ax1.plot(x_vector, decay(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2]), 'r', linewidth=linewidth)
    else:
        ax1.plot(x_vector, ramsey(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2], kwargs['fitted_pars'][3], kwargs['fitted_pars'][4]), 'r', linewidth=linewidth)
    ax1.set_title('Ramsey %03d' % (kwargs['iteration']))
    textstr = format_textstr_ramsey(kwargs)
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)
    plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=8)


def plot_echo(x_vector, y_vector, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_vector, y_vector, '-o', markersize=3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Separation ($\mu$s)')
    ax.plot(x_vector, decay(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2]), 'r')
    textstr = format_textstr_echo(kwargs)
    ax.set_title('Echo Measurement %03d' % (kwargs['iteration']))
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)
    plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=8)


def plot_cavity_reset(x_vector, y_vector, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_vector, y_vector, '-o', markersize=3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Separation ($\mu$s)')
    ax.plot(x_vector, decay(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2]), 'r')
    textstr = format_textstr_cavity_reset(kwargs)
    ax.set_title('Echo Measurement %03d' % (kwargs['iteration']))
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)
    plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=8)


def plot_cavity_cooling(x_vector, y_vector, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_vector, y_vector, '-o', markersize=3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Separation ($\mu$s)')
    ax.plot(x_vector, decay(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2]), 'r')
    textstr = format_textstr_cavity_cooling(kwargs)
    ax.set_title('Echo Measurement %03d' % (kwargs['iteration']))
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)
    plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=8)


def plot_ringdown(x_vector, y_vector, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_vector, y_vector, '-o', markersize=3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Delay ($\mu$s)')
    ax.plot(x_vector, decay(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2]), 'r')
    textstr = format_textstr_ringdown(kwargs)
    ax.set_title('resonator ring down %03d' % (kwargs['iteration']))
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)
    plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=8)


def plot_T1(x_vector, y_vector, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_vector, y_vector, '-o', markersize=3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Delay ($\mu$s)')
    ax.plot(x_vector, decay(x_vector, kwargs['fitted_pars'][0], kwargs['fitted_pars'][1], kwargs['fitted_pars'][2]), 'r')
    textstr = format_textstr_T1(kwargs)
    ax.set_title('T1 Measurement %03d' % (kwargs['iteration']))
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)
    plt.tick_params(axis='both', direction='in', bottom=True, top=True, left=True, right=True, size=8)

def format_textstr_p_rabi(kwargs):
    return (
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.6f} GHz\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
        f'rr atten={kwargs["rr_atten"]:.1f} db\n'
        f'gauss len={kwargs["gauss_len"]} ns\n'
        f'$A_{{\pi/2}}$={kwargs["fitted_pars"][1] / 4 * kwargs["gauss_amp"]:.3f} V'
    )

def format_textstr_rabi(kwargs):
    return (
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.4f} GHz\n'
        '$T_{X90}$' +f'={round(kwargs["fitted_pars"][1] / 4, 1):.1f} ns\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
        f'rr atten={kwargs["rr_atten"]:.1f} db\n'
        f'contrast={2*kwargs["fitted_pars"][0]:.2f} mV\n'
        f'gauss amp={kwargs["gauss_amp"]:.4f}'
    )

def format_textstr_ramsey(kwargs):
    return (
        '$T_{X90}$' +f'={kwargs["X90_len"]:.1f} ns\n'
        '$A_{X90}$'+f'={kwargs["X90_amp"]:.3f} V\n'
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.6f} GHz\n'
        f'$\Delta$ = {kwargs["fitted_pars"][1]:.2f} MHz\n'
        f'$T_2$ = {kwargs["fitted_pars"][3]:.3f} $\mu$s\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
        f'rr atten={kwargs["rr_atten"]:.1f} db'
    )

def format_textstr_echo(kwargs):
    return (
        '$T_{X90}$' +f'={kwargs["X90_len"]:.1f} ns\n'
        '$A_{X90}$'+f'={kwargs["X90_amp"]:.3f} V\n'
        '$A_{X180}$'+f'={kwargs["X180_amp"]:.3f} V\n'
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.4f} GHz\n'
        f'$T_2$={kwargs["fitted_pars"][1]:.2f}$\mu$s\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
    )

def format_textstr_cavity_reset(kwargs):
    return (
        f'$T_{{\pi/2}}$={kwargs["pi2Width"]:.1f} ns\n'
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.4f} GHz\n'
        f'$T_2$={kwargs["fitted_pars"][1]:.2f}$\mu$s\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
        f'$T_2error$ = {kwargs["error"][1]:.2f} us\n'
        f'$flux$= {kwargs["flux"]:.3f} mA\n'
        f'$amp ffl scale$={kwargs["amp_ffl_scale"]:.2f}\n'
        f'$amp rr scale$={kwargs["amp"]:.2f}\n'
        f'$ffl len$={kwargs["ffl_len"]}\n'
        f'$\omega_{{ffl}}$={kwargs["fflDriveFreq"] * 1e-9:.3f} GHz\n'
        f'$(ffl,rr) atten$=({kwargs["ffl_atten"]:.1f}, {kwargs["rr_atten"]:.1f}) db'
    )

def format_textstr_cavity_cooling(kwargs):
    return (
        f'$T_{{\pi/2}}$={kwargs["pi2Width"]:.1f} ns\n'
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.4f} GHz\n'
        f'$A_d$ = {kwargs["qb_power"]:.2f} V\n'
        f'$T_2$={kwargs["fitted_pars"][1]:.2f}$\mu$s\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
        f'$T_2error$ = {kwargs["error"][1]:.2f} us\n'
        f'$flux$= {kwargs["flux"]:.3f} mA\n'
        f'$amp ffl scale$={kwargs["amp_ffl_scale"]:.2f}\n'
        f'$amp rr scale$={kwargs["amp"]:.2f}\n'
        f'$\omega_{{ffl}}$={kwargs["fflDriveFreq"] * 1e-9:.3f} GHz\n'
        f'$(ffl,rr) atten$=({kwargs["ffl_atten"]:.1f}, {kwargs["rr_atten"]:.1f}) db'
    )

def format_textstr_ringdown(kwargs):
    return (
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.4f} GHz\n'
        f'$T_{{ringdown}}$ = {kwargs["fitted_pars"][1]:.4f} $\mu$s\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
        f'$flux$= {kwargs["flux"]:.3f} mA\n'
        f'$amp ffl scale$={kwargs["amp_ffl_scale"]:.2f}\n'
        f'$T_{{ring}}error$ = ({kwargs["error"][1]:.2f})us\n'
        f'$(ffl,rr) atten$=({kwargs["ffl_atten"]:.1f}, {kwargs["rr_atten"]:.1f}) db'
    )

def format_textstr_T1(kwargs):
    return (
        f'$T_{{\pi/2}}$={kwargs["X180_len"]:.1f} ns\n'
        f'$A_{{\pi}}$={kwargs["X180_amp"]:.3f}\n'
        f'$\omega_d$ = {kwargs["qubitDriveFreq"] * 1e-9:.4f} GHz\n'
        f'$T_1$ = {kwargs["fitted_pars"][1]:.3f} $\mu$s\n'
        '$\hat{n}$'+ f'= {kwargs["nAverages"]}\n'
        f'rr atten= {kwargs["rr_atten"]:.1f} db'
    )


def rabi(x, amp,period,phase,offset):
    return amp*np.cos(2*pi*x/period+phase)+offset

def ramsey(x,amp,f,phase,tau,offset):
    return amp*np.cos(2*pi*f*x+phase)*np.exp(-x/tau)+offset

def beats(x,amp,f1,f2,phase1,phase2,tau,offset):
    return amp*np.cos(pi*(f1+f2)*x+phase1)*np.cos(pi*(f2-f1)*x+phase2)*np.exp(-x/tau)+offset

def decay(x,amp,tau,offset):
    return amp*np.exp(-x/tau)+offset

def mod_cos(x,amp,B0,nu,phi1,phi2,tau,offset):
    return amp*np.cos(B0/nu*np.sin(2*np.pi*nu*x+phi1)+phi2)*np.exp(-x/tau)+offset
    # return amp*np.cos(np.cos(nu*x)*f*x)*np.exp(-x/tau)+offset

def mod_dec(x,amp1,f1,phi1,tau1,amp2,phi2,tau2,offset):
    return amp1*np.cos(2*np.pi*f1*x+phi1)*np.exp(-x/tau1)+ amp2*np.sin(2*np.pi*f1*x+phi2)*np.exp(-x/tau2)+offset
    # return amp*np.cos(np.cos(nu*x)*f*x)*np.exp(-x/tau)+offset

def extract_freq(t_vector, y_vector, dt, plot=0):
    """
    Extracts the dominant frequency from a time-domain signal using FFT.

    Parameters:
    t_vector (array-like): Time vector.
    y_vector (array-like): Signal vector.
    dt (float): Time step between samples.
    plot (int, optional): If 1, plots the power spectral density. Default is 1.

    Returns:
    float: Dominant frequency in the signal.
    """
    N = len(t_vector)
    yf = scy.fft.fft(y_vector - np.mean(y_vector))
    xf = scy.fft.fftfreq(N, dt)[:N // 2]
    psd = 2.0 / N * np.abs(yf[:N // 2])
    index_max = np.argmax(psd)

    if plot == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xf, psd)
        ax.set_xlabel('Frequency (MHz)')
        ax.set_ylabel('Power')

    return xf[index_max]

def get_envelope(sig,dt, distance):
    # split signal into negative and positive parts
    sig = sig - np.mean(sig)
    u_x = np.where(sig > 0)[0]
    l_x = np.where(sig < 0)[0]
    u_y = sig.copy()
    u_y[l_x] = 0
    l_y = -sig.copy()
    l_y[u_x] = 0

    # find upper and lower peaks
    u_peaks, _ = scipy.signal.find_peaks(u_y, distance=distance)
    l_peaks, _ = scipy.signal.find_peaks(l_y, distance=distance)
    # use peaks and peak values to make envelope
    u_x = u_peaks
    u_y = sig[u_peaks]
    l_x = l_peaks
    l_y = sig[l_peaks]

    # add start and end of signal to allow proper indexing
    end = len(sig)
    u_x = np.concatenate(([0],u_x, [end]))*dt
    u_y = np.concatenate(([sig[0]],u_y, [sig[-1]]))
    l_x = np.concatenate(([0],l_x, [end]))*dt
    l_y = np.concatenate(([min(sig)],l_y, [np.mean(sig)]))
    # create envelope functions
    u = scipy.interpolate.interp1d(u_x, u_y,kind='cubic',fill_value="extrapolate")
    # l = scipy.interpolate.interp1d(l_x, l_y,kind='cubic')
    return u

def get_envelope_LPF(x,sig):

    N = len(sig)
    Tmax = x[-1]
    cutoff = 100e6
    fs = N/Tmax

    env = butter_lowpass_filter(sig, cutoff, fs)
    plt.plot(x,env)
    plt.show()

    return env

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    sos = butter(order, normal_cutoff, btype='low', output = 'sos', analog=False)
    return sos

def butter_lowpass_filter(data,cutoff,fs):
    sos = butter_lowpass(cutoff, fs, order=5)
    y = sp.signal.sosfilt(sos, data)
    return y

def Volt2dBm(data):

    return 10*np.log10(1e3*data**2/50)

def Watt2dBm(x):
    '''
    converts from units of Watts to dBm
    '''
    return 10.*np.log10(x*1000.)


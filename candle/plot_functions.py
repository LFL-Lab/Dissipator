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

#%% spec_plot
def spec_plot(freq,I,Q,attenuation=-30,df=0.1e6,plot='mag',element='resonator',fwhm=0,fc=0,iteration=1,find_peaks=False):

    freq = freq*1e-9
    I = I*1e3
    Q = Q*1e3
    mag = np.abs(I+1j*Q)

    phase = np.unwrap(np.angle(I+1j*Q))
    if element == 'qubit' and find_peaks:
        sigma = np.std(mag)
        peaks,_ = scy.signal.find_peaks(mag,height=np.mean(mag)+5*sigma,distance=200,width=10)
        print(peaks)
        try:
            print(f'Peaks at: {round(freq[peaks[0]],5)} GHz, {round(freq[peaks[1]],5)}')
        except:
            print('Peaks not found or do not exist.')

    fig = plt.figure(figsize=(8,8))

    if element == 'qubit':

        # I data
        ax1 = fig.add_subplot(221)
        ax1.plot(freq,I,'-o', markersize = 3, c='C0')
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('I (mV)')
        # Q data
        ax1 = fig.add_subplot(222)
        ax1.plot(freq,Q,'-o', markersize = 3, c='C0')
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('Q (mV)')
        # phase data
        ax1 = fig.add_subplot(223)
        ax1.plot(freq,phase,'-o', markersize = 3, c='C0')
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('Phase (rad)')
        # Power data
        ax1 = fig.add_subplot(224)
        ax1.plot(freq,mag,'-o', markersize = 3, c='C0')
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('Magnitude (mV)')

    elif element == 'resonator':
        # Power data
        ax1 = fig.add_subplot(211)
        ax1.plot(freq,mag,'-o', markersize = 3, c='C0')
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('Magnitude (mV)')
        ax2 = fig.add_subplot(212)
        ax2.plot(freq,phase,'-o', markersize = 3, c='C0')
        ax2.set_xlabel('Frequency (GHz)')
        ax2.set_ylabel('Phase (deg)')

    if element == 'resonator':
        txt = f'$\omega_c$ = {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-6:.3f} MHz\n$\kappa$ = {2*np.pi*fwhm*1e-6:.3f} MHz\nReadout attenuation: {attenuation} dB\ndf = {df*1e-3} kHz'
    elif element == 'qubit':
        if len(peaks) == 2:
            txt = f'$\omega_{10}$ = {round(freq[peaks[1]],3)} GHz\n$\omega_{12}$ = {round(freq[peaks[0]],3)} GHz\n$\\alpha$ = {round((freq[peaks[0]]-freq[peaks[1]])*1e3,1)} MHz'
        else:
            txt = ''
    plt.gcf().text(1, 0.15, txt, fontsize=14)
    plt.title(f'{element} spectroscopy {iteration}')
    plt.tight_layout()

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
    fig,ax = plt.subplots(figsize=(4,3), dpi=300)
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
    plt.tight_layout()

    return hm;

#%% plot_single_shot
def plot_single_shot(datadict,axes=0):

    datadict = {key: value*1e3 for key,value in datadict.items()} # convert to mV
    datadict = {key: value.tolist() for key,value in datadict.items()} # convert to list

    states = []
    for key,value in datadict.items():
        print(key+':'+str(len(value))+'\n')
    [states.append(r'$|g\rangle$') for i in range(len(datadict['I']))]
    [states.append(r'$|e\rangle$') for i in range(len(datadict['Iexc']))]
    data = {
            'I [mV]':   np.hstack((datadict['I'],datadict['Iexc'])),
            'Q [mV]':   np.hstack((datadict['Q'],datadict['Qexc'])),
            'States':   states
                }
    dataF = pd.DataFrame(data=data)
    plot = sns.jointplot(data=dataF, x='I [mV]',y='Q [mV]',hue='States',ax=axes,space=0)

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


#%% fit_res
def fit_res(f_data,z_data,res_type='notch'):
    fc = f_data[np.argmin(z_data)]
    if res_type == 'notch':
        z_data = -z_data-min(-z_data)
        idx = np.argwhere(np.diff(np.sign(z_data - 0.5*max(z_data)))).flatten()
        fwhm = f_data[idx[1]] - f_data[idx[0]]

    return fc,fwhm

#%% fit_data
def fit_data(x_vector,y_vector,sequence='rabi',dt=0.01,fitFunc='',verbose=0):

    '''
    fit experimental data

    sequence:          'Rabi','ramsey', 'T1' or 'T2'
    x_vector:           time data
    y_vector:           voltage data
    dt:                 sequence stepsize. Used for extracting the frequency of the data
    '''
    x_vector = x_vector
    y_vector = y_vector*1e3

    amp = (max(y_vector)-min(y_vector))/2
    offset = np.mean(y_vector)

    if sequence == "rabi":
        fitFunction = rabi
        period = 1e3/(extract_freq(x_vector*1e3, y_vector, dt,plot=0))
        print('Period Initial Guess: %.1f ns'%(period))
        phase = pi
        x_vector = x_vector*1e3
        lb = [0.1*amp,0.1*period,0,-2*abs(offset)]
        ub = [10*amp,10*period,2*pi,2*abs(offset)]
        p0 = [amp,period,phase,offset]

    elif sequence == "ramsey":
        f = extract_freq(x_vector, y_vector,dt,plot=0)
        # print('Initial Guess for Freq:%.4f MHz'%(f))
        if x_vector[-1] > 20:
            tau = 30
        else:
            tau = 2
        phi = 0
        amp = abs(amp)
        # try:
        if fitFunc != 'envelope':
            p0 = [amp,f,phi,tau,offset]
            lb = [0.75*amp,0.1*f,-pi,0.01,-2*abs(offset)]
            ub = [2*amp,2*f,pi,100,2*abs(offset)]
            fitFunction = ramsey
            # fitted_pars, covar = scy.optimize.curve_fit(fitFunction, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)
        elif fitFunc == 'envelope':
            tau = 1
            env = get_envelope(y_vector, dt, distance=100)
            env = env(x_vector) + offset
            # env = get_envelope_LPF(x_vector, y_vector)*1e-3
            p0 = [amp,tau,offset]
            if offset < 0:
                p0 = [amp,tau,offset]
                lb = [0.95*amp,0.1,2*offset]
                ub = [2*amp,15,0.5*offset]
            elif offset >= 0:
                p0 = [amp,tau,offset]
                lb = [0.9*amp,0.1,0.9*offset]
                ub = [1.1*amp,15,1.1*offset]
            fitFunction = decay
            y_vector = env
            # fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, env,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)
    elif sequence == "echo":
        if x_vector[-1] < 10:
            tau = 2
            tau_ub = 20
        else:
            tau = 20
            tau_ub = 300
        amp = y_vector[0] - y_vector[-1]
        p0 = [amp,tau,offset]
        if offset < 0:
            lb = [0.95*amp,0.1,1.05*offset]
            ub = [1.05*amp,tau_ub,0.95*offset]
        elif offset >= 0:
            lb = [0.95*amp,0.1,0.95*offset]
            ub = [1.05*amp,tau_ub,1.05*offset]
        fitFunction = decay
        # fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=6000)
    elif sequence == "T1":
        tau = 2
        amp = y_vector[0] - y_vector[-1]
        if amp < 0:
            p0 = [amp,tau,offset]
            lb = [10*amp,0.1,-2*abs(offset)]
            ub = [0.5*amp,300,2*abs(offset)]
        elif amp >= 0:
            p0 = [amp,tau,offset]
            lb = [0.5*amp,0.1,-2*abs(offset)]
            ub = [10*amp,300,2*abs(offset)]
        fitFunction = decay
        # fitted_pars, covar = scy.optimize.curve_fit(, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=6000)

    fitted_pars, covar = scy.optimize.curve_fit(fitFunction, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=40e3)
    error = np.sqrt(abs(np.diag(covar)))

    if verbose == 1:
        print('-'*100)
        print('Lower Bounds:',np.around(lb,1))
        print('Initial Guess:',np.around(p0,1))
        print('Upper Bounds:',np.around(ub,1))
        print('Best Fit Pars:',np.around(fitted_pars,1))
        print('Error:',np.around(error,1))
        print('-'*100)
    else:
        pass

    return fitted_pars,error


#%% plot_data
def plot_data(x_vector,y_vector,sequence='rabi',qubitDriveFreq=3.8e9,amplitude_hd=1,
                              pi2Width='',nAverages=1,
                              integration_length=2e-6,cav_resp_time=5e-6,stepSize=5e-6, iteration = 1,
                              Tmax=5e-6,measPeriod=5e-6,active_reset=False,
                              fitted_pars=np.zeros(7),plot_mode=0,rr_IF=5e6,fitFunc='',savefig=True):

    x_vector = x_vector*1e3
    y_vector = y_vector*1e3

    if sequence == "a-rabi":
        fig, ax = plt.subplots()
        ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Pulse Amplitude (mV)')
        # ax.plot(x_vector*1e3,rabi(x_vector*1e3, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
        ax.set_title('Power Rabi Measurement %03d'%(iteration))
        # textstr = '$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_{\pi/2}$ = %.1f ns\n$\hatn$ = %d'%(qubitDriveFreq*1e-9,amplitude_hd,round(fitted_pars[1]/4,1),nAverages)
        # plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    if sequence == "t-rabi":
        fig, ax = plt.subplots()
        ax.plot(x_vector*1e3, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Pulse Duration (ns)')
        ax.plot(x_vector*1e3,rabi(x_vector*1e3, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
        ax.set_title('Rabi Measurement %03d'%(iteration))
        textstr = '$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_{\pi/2}$ = %.1f ns\n$\hatn$ = %d'%(qubitDriveFreq*1e-9,amplitude_hd,round(fitted_pars[1]/4,1),nAverages)
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    elif sequence == "ramsey":

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        fontSize = 16
        tickSize = 11
        markersize= 4
        linewidth = 2
        ax1.plot(x_vector, y_vector, 'o', markersize = markersize, c='C0')
        ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
        ax1.set_xlabel('Pulse Separation ($\mu$s)')
        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(tickSize)
        if fitFunc == 'envelope':
            ax1.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r',linewidth=linewidth)
        else:
            # print(len(fitted_pars))
            ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
        # ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
        # textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.1f mV\n$\\tau_k$ = %.3f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],mu*1e3,AC_freq*1e-9,sigma*1e3,B0*1e3,nu,nAverages)
        # plt.gcf().text(0.95, 0.15, textstr, fontsize=fontSize)
        # ax1.set_title('Ramsey Measurement %03d' %(iteration),size=fontSize)

    elif sequence == "echo":

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Pulse Separation ($\mu$s)')
        ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
        textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_2$=%.2f$\mu$s\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],nAverages)
        ax.set_title('Echo Measurement %03d' %(iteration))
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)


    elif sequence == "T1":
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Delay ($\mu$s)')
        ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
        textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_1$=%.2f $\mu$s\n$\mu$ = %.3f V\n\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],nAverages)
        ax.set_title('T1 Measurement %03d' %(iteration))
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    if savefig:
        plt.savefig(f'D:\\weak_measurements\\{sequence}\\fig_{iteration:03d}.png',dpi='figure')

    return fig

def extract_data(sweep,B0,nu,tauk,meas_device='CandleQubit_6',nMeasurements=100,nBackMeasurements=100,fileformat='new',sequence='ramsey',iteration=1):

    # get data
    if fileformat == 'old':
        filename = 'B0_%d_uV_nu_%d_kHz_tau_%d_ns' %(round(B0*1e6),round(nu*1e3),round(tauk*1e3))
        datafile = "E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep,filename)
        tdata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=3,header=None,nrows=1).dropna(axis='columns').to_numpy(np.float64)[0]
        tdata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=5,header=None,nrows=1).dropna(axis='columns').to_numpy(np.float64)[0]
        ydata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7,header=None,nrows=nBackMeasurements).dropna(axis='columns').to_numpy(np.float64)
        ydata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7+2*(nBackMeasurements+1),header=None,nrows=nMeasurements).dropna(axis='columns').to_numpy(np.float64)
        # get exp parameters
        pars = pd.read_csv(datafile,on_bad_lines='skip',header=[0],nrows=1)
        keys = pars.keys()
        values = pars.values
        dictionary = {"mu":loads(pars.loc[0].at['AC_pars'])[0],
                      "sigma":loads(pars.loc[0].at['AC_pars'])[1],
                      "B0":loads(pars.loc[0].at['RT_pars'])[0],
                      # "nu": loads(pars.loc[0].at['RT_pars'])[1],
                      "tauk":loads(pars.loc[0].at['RT_pars'])[1]}
    elif fileformat == 'new':
        # cols_background = [i for i in range(0,94)]
        # cols = [i for i in range(0,250)]
        filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns' %(round(B0*1e6),round(nu*1e3),round(tauk*1e3))
        # filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns_%d' %(round(B0*1e6),round(nu*1e3),round(tauk*1e3),iteration)
        datafile = "E:\\generalized-markovian-noise\\%s\\sweep_data\\%s\\%s\\data\\data_%s.csv"%(meas_device,sequence,sweep,filename)
        tdata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=3,header=None,nrows=1).to_numpy(np.float64)[0]
        tdata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=5,header=None,nrows=1).to_numpy(np.float64)[0]
        ydata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7,header=None,nrows=nBackMeasurements).to_numpy(np.float64)
        ydata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7+nBackMeasurements+1,header=None,nrows=nMeasurements).to_numpy(np.float64)
        # get exp parameters
        pars = pd.read_csv(datafile,on_bad_lines='skip',header=[0],nrows=1)
        keys = pars.keys()
        values = pars.values
        dictionary = dict(zip(keys,values[0]))
    #print(dictionary)

    return tdata_background, ydata_background, tdata, ydata, dictionary


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

def extract_freq(t_vector,y_vector,dt,plot=0):
    N = len(t_vector)
    dt = dt*1e6
    yf = scy.fft.fft(y_vector-np.mean(y_vector))
    xf = scy.fft.fftfreq(N,dt)[:round(N/2)]
    # print(len(xf))
    psd = 2.0/N * np.abs(yf[:round(N/2)])
    # print(len(psd))
    # print(psd)
    index_max = np.argmax(psd)
    if plot == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xf,psd)
        ax.set_xlabel('Frequency (MHz)')
        ax.set_ylabel('Power')
    # print(index_max)

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


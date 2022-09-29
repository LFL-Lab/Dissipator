# -*- coding: utf-8 -*-
"""
Created on Mon May 24 13:13:39 2021

@author: lfl
"""
import matplotlib.pyplot as plt


import numpy as np
import pandas as pd
import scipy as sp
import scipy as scy
# from matplotlib import cm
import sympy as sy
import csv
import itertools
from scipy.interpolate import interp1d
import scipy.fftpack
import time
import os
from json import loads
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LightSource
from types import SimpleNamespace
pi=np.pi
import zhinst.utils as ziut
import seaborn as sns; sns.set() # styling
from matplotlib.ticker import FormatStrFormatter
from experiment_funcs import gen_tel_noise
from scipy.signal import butter,lfilter,freqz
# import imageio
from experiment_funcs import convertB0_mV_to_kHz,convertB0_kHz_to_mV
sns.set_style('ticks')
# plt.rcParams['figure.dpi'] = 300

def spec_plot(freq,I,Q,readout_power=-30,qubit_drive_amp=0.2):

    freq = freq*1e9
    I = I*1e3
    Q = Q*1e3
    mag = np.abs(I*I.conjugate()+Q*Q.conjugate())

    phase = np.unwrap(np.angle(I+1j*Q))
    # sigma = np.std(I)
    # peak = scy.signal.find_peaks(I,height=np.mean(I)+sigma,distance=100)
    # # print(peak[0])
    # peaks = peak[0]

    fig = plt.figure(figsize=(10,7))

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
    ax1.set_ylabel('Mag (mV)')

    # textstr = "$P_r = %.1f$ dBm\n Qubit Wfm Amp = %.1f mV" %(readout_power,qubit_drive_amp*1e3)
    # plt.gcf().text(1, 0.25, textstr, fontsize=14)

    # txt = "\n"
    # for i in peaks:
    #     txt += "%.4f GHz \n" %(freq[i]*1e-9)
    # print('Peaks are at: %s'%txt)

def fit_data(x_vector,y_vector,sequence='rabi',dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,mu=0,sigma=0,B0=0,nu=0,tauk=0,
                 pi2Width=50e-9,piWidth_Y=0,fitting = 1, plot=1,save_fig = 0,iteration=1,nAverages=1,sampling_rate=2.4e9,
                 integration_length=2e-6,AC_freq=5e9,cav_resp_time=5e-6,sweep=0,stepSize=5e-6,Tmax=5e-6,measPeriod=5e-6,nSteps=101,
                 prePulseLength=1500e-9,postPulseLength=1500e-9,threshold=10e-3,active_reset=False,rr_IF=5e6,source=0,noiseType=1,fitFunc='',
                 pipulse_position='',phi=0,axis='X',wk=0,noise_rate=1,verbose=1,simulation=0):

    '''
    fit experiment data

    sequence:          'Rabi', 'T1' or 'T2'
    complex_amplitude:  complex amplitude from the measurement
    x_vector:           x data
    fitting:     0:     do not fit
                  1:     do fit
    save_fig:    0:     do not save
                  1:     do save
    '''
    if simulation == 0:
        x_vector = x_vector*1e6
        y_vector = y_vector*1e3
    elif simulation == 1:
        pass

    amp = (max(y_vector)-min(y_vector))/2
    offset = np.mean(y_vector)

    if sequence == "rabi":
        fitFunction = rabi
        period = 1e3/(extract_freq(x_vector*1e3, y_vector, dt,plot=0))
        # print('Period Initial Guess: %.1f ns'%(period))
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
        if noiseType == 1 and fitFunc != 'envelope':
            p0 = [amp,f,phi,tau,offset]
            lb = [0.75*amp,0.1*f,-pi,0.01,-2*abs(offset)]
            ub = [2*amp,2*f,pi,100,2*abs(offset)]
            fitFunction = ramsey
            # fitted_pars, covar = scy.optimize.curve_fit(fitFunction, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)
        elif noiseType == 2 and fitFunc == 'modDampedCos':
            tau = 5
            B0 = 4.02e-2*B0*1e3+1.53e-3
            p0 = [amp,B0,phi,phi,tau,offset]
            lb = [0.96*amp,0.1*B0,0,0,0.1,-2*abs(offset)]
            ub = [1.05*amp,10*B0,2*pi,2*pi,100,2*abs(offset)]
            fitFunction = lambda x,amp,B0,phi1,phi2,tau,offset: mod_cos(x,amp,B0,nu*1e-3,phi1,phi2,tau,offset)
            # fitted_pars, covar = scy.optimize.curve_fit(lambda x,amp,B0,phi,tau,offset: mod_cos(x,amp,B0,nu*1e-3,phi,tau,offset), x_vector, y_vector, p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)

        elif noiseType == 2 and fitFunc == 'modDecayPlusCos':
            B0 = 4.02e-2*B0*1e3+1.53e-3
            amp1 = amp
            amp2 = amp
            phi2 = pi/2
            tau = 5
            tau2 = 5
            p0 = [amp1,B0,phi,tau,amp2,phi2,tau2,offset]
            if offset < 0:
                offset_lb = 2*offset
                offset_ub = 0.5*offset
            elif offset > 0:
                offset_lb = 0.5*offset
                offset_ub = 2*offset
            lb = [0.75*amp, 0, -pi, 0, 0.5*amp, 0, 0.001, offset_lb]
            ub = [1.5*amp, 100*B0, pi, 100, 2*amp,  2*pi, 100, offset_ub]
            fitFunction = mod_dec
            # fitted_pars, covar = scy.optimize.curve_fit(mod_dec, x_vector, y_vector, p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)

        elif (noiseType == 2 or noiseType == 1) and fitFunc == 'envelope':
            tau = 1
            if simulation == 0:
                env = get_envelope(y_vector*1e-3, dt*1e6, distance=80)
                env = 1e3*env(x_vector) + offset
            elif simulation == 1:
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
        # p0 = [amp,tau,offset]
        # lb = [0.75*amp,0.1,2*offset]
        # ub = [2*amp,200,0.5*offset]
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

def plot_data(awg,x_vector,y_vector,error=0,sequence='rabi',dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,mu=0,sigma=0,B0=0,nu=0,tauk=0,
                 pi2Width=50e-9,piWidth_Y=0,fitting = 1, plot=1,save_fig = 0,iteration=1,nAverages=1,sampling_rate=2.4e9,
                 integration_length=2e-6,AC_freq=5e9,cav_resp_time=5e-6,sweep=0,stepSize=5e-6,Tmax=5e-6,measPeriod=5e-6,nSteps=101,
                 prePulseLength=1500e-9,postPulseLength=1500e-9,threshold=10e-3,active_reset=False,fitted_pars=np.zeros(7),pi_pulse=50,
                 fit_single_par_point=0,plot_mode=0,rr_IF=5e6,source=0,noiseType=1,fitFunc='',pipulse_position='',phi=0,axis='X',wk=0,noise_rate=1,
                 simulation=0):

    if fit_single_par_point == 0 and simulation == 0:
        x_vector = x_vector*1e6
        y_vector = y_vector*1e3
    elif fit_single_par_point == 1 and simulation == 0:
        x_vector[0] = x_vector[0]*1e6
        y_vector[0] = y_vector[0]*1e3
        x_vector[1] = x_vector[1]*1e6
        y_vector[1] = y_vector[1]*1e3
    elif simulation == 1:
        pass

    if sequence == "rabi":
        fig, ax = plt.subplots()
        ax.plot(x_vector*1e3, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Pulse Duration (ns)')
        ax.plot(x_vector*1e3,rabi(x_vector*1e3, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
        ax.set_title('Rabi Measurement %03d'%(iteration))
        textstr = '$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_{\pi/2}$ = %.1f ns\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\hatn$ = %d'%(qubitDriveFreq*1e-9,amplitude_hd,round(fitted_pars[1]/4,1),mu*1e3,AC_freq*1e-9,nAverages)
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    elif sequence == "ramsey":

        if simulation == 1:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            fontSize = 16
            tickSize = 11
            markersize= 4
            linewidth = 2
            ax1.plot(x_vector, y_vector, 'bo', markersize = markersize,label='Data')
            ax1.set_ylabel('Fidelity',fontsize=fontSize)
            ax1.set_xlabel('Time ($\mu$s)')

            for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                label.set_fontsize(tickSize)
            if fitFunc == 'envelope':
                ax1.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r',linewidth=linewidth,label='Fit')
            else:
                # print(len(fitted_pars))
                ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth,label='Fit')
            ax1.legend()
            ax1.set_title('$B_0$ = %d kHz | $\\nu$ = %d kHz | $\\tau_k$ = %1.f $\mu$s'%(B0,nu,tauk))

        if (B0 != 0 or mu != 0 or sigma != 0) and simulation == 0: # if noise is injected, sz or sx, plot data and noise instances
            if fit_single_par_point == 1: # for plotting data averaged over many noise realizations post par sweep
                fig = plt.figure(figsize=((14,8)))
                ax1 = fig.add_subplot(111)
                # Plot data
                fontSize = 22
                tickSize = 20
                markersize= 10
                linewidth = 6

                B0 = int(convertB0_mV_to_kHz(B0*1e3))
                # plot data
                ax1.plot(x_vector[0], y_vector[0], 'o', markersize = markersize, c='C0',label='$B_0 = 0$')
                if nu != 0:
                    if tauk > 100:
                        label='$B_0$ = %d kHz | $\\nu$ = %.1f kHz | $\\tau_k/\\tau_0$ = $\infty$'%(B0,nu)
                    else:
                        label='$B_0$ = %d kHz | $\\nu$ = %.1f kHz | $\\tau_k/\\tau_0$ = %.1f'%(B0,nu,tauk/fitted_pars[0][1])
                elif nu == 0:
                    if tauk > 100:
                        label='$B_0$ = %d kHz | $\\tau_k/\\tau_0$ = $\infty$'%(B0)
                    else:
                        label='$B_0$ = %d kHz | $\\tau_k/\\tau_0 = %.1f$'%(B0,tauk/fitted_pars[0][1])
                ax1.plot(x_vector[1], y_vector[1], 'o', markersize = markersize, c='k',label=label)

                # plot fits
                # background fit
                ax1.plot(x_vector[0],decay(x_vector[0], fitted_pars[0][0], fitted_pars[0][1], fitted_pars[0][2]),'r',linewidth=linewidth,label='$\\tau_0 = %.1f \mu s$'%(fitted_pars[0][1]))

                if fitted_pars[1][1] < 0.01:
                    fitted_pars[1][1] = 0


                # background + generalized markovian data fit
                if noiseType == 1 and fitFunc != 'envelope':
                    # err = fitted_pars[1][3]/fitted_pars[0][3]*np.sqrt((error[0]/fitted_pars[0][3])**2+(error[1]/fitted_pars[1][3])**2)
                    ax1.plot(x_vector[1],ramsey(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], fitted_pars[1][2],fitted_pars[1][3],fitted_pars[1][4]),'g',linewidth=linewidth,label='$\omega = 2\pi\\times%.1f$ MHz, $\\tau/\\tau_0$ = %.2f '%(fitted_pars[1][1],fitted_pars[1][3]/fitted_pars[0][1]))
                elif noiseType == 1 and fitFunc == 'envelope':
                    # print(fitted_pars)
                    # err = fitted_pars[1][3]/fitted_pars[0][3]*np.sqrt((error[0]/fitted_pars[0][3])**2+(error[1]/fitted_pars[1][3])**2)
                    ax1.plot(x_vector[1],decay(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], fitted_pars[1][2]),'g',linewidth=linewidth,label='$\omega = 2\pi\\times%.1f$ MHz, $\\tau/\\tau_0$ = %.2f '%(fitted_pars[1][1],fitted_pars[1][1]/fitted_pars[0][1]))
                elif noiseType == 2 and fitFunc == 'modDampedCos':
                    ax1.plot(x_vector[1], mod_cos(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], nu*1e-3, fitted_pars[1][2], fitted_pars[1][3],fitted_pars[1][4],fitted_pars[1][5]),'g',linewidth=linewidth,label='$\\tau/\\tau_0$ = %.1f'%(fitted_pars[1][4]/fitted_pars[0][3]))
                elif noiseType == 2 and fitFunc == 'modDecayPlusCos':
                    ax1.plot(x_vector[1], mod_dec(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], fitted_pars[1][2], fitted_pars[1][3],fitted_pars[1][4],fitted_pars[1][5],fitted_pars[1][6],fitted_pars[1][7]),'g',linewidth=linewidth,label='$\\tau_1/\\tau_0$ = %.1f, $\\tau_2/\\tau_0$ = %.1f, $f_1$ = %.1f MHz'%(fitted_pars[1][3]/fitted_pars[0][1],fitted_pars[1][6]/fitted_pars[0][3],fitted_pars[1][1]))
                elif noiseType == 2 and fitFunc == 'envelope':
                    ax1.plot(x_vector[1],decay(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], fitted_pars[1][2]),'g',linewidth=linewidth,label='$\omega = 2\pi\\times%.1f$ MHz, $\\tau/\\tau_0$ = %.1f'%(B0*1e-3,fitted_pars[1][1]/fitted_pars[0][1]))

                ax1.legend(loc='upper right',prop={'size':22})
                ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
                ax1.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
                for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                    label.set_fontsize(tickSize)

                if x_vector[0][-1] < x_vector[1][-1]:
                    ax1.set_xlim([0,x_vector[0][-1]])
                else:
                    ax1.set_xlim([0,x_vector[1][-1]])
                # textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.5f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$\\tau$/$\\tau_0$=%.1f\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[0][1],fitted_pars[1][3]/fitted_pars[0][3],mu*1e3,AC_freq*1e-9,sigma*1e3,nAverages)
                # plt.gcf().text(0.95, 0.15, textstr, fontsize=fontSize)
            elif fit_single_par_point == 0: #for plotting non-par sweep data and par sweep data during the sweep
                fig = plt.figure(figsize=(20,18))
                ax1 = fig.add_subplot(2,2,(1,2))
                # Plot data
                fontSize = 36
                tickSize = 24
                markersize= 12
                linewidth = 6
                ax1.plot(x_vector, y_vector, 'o', markersize = markersize, c='C0')
                ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
                for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                    label.set_fontsize(tickSize)
                # Retrieve and plot telegraph noise instrance
                waveforms = ziut.parse_awg_waveform(awg.get('/dev8233/awgs/0/waveform/waves/0')['dev8233']['awgs']['0']['waveform']['waves']['0'][0]['vector'],channels=2)
                t_arr = np.linspace(0,x_vector[-1],len(waveforms[0]))
                ax2 = fig.add_subplot(2,2,3)
                ax2.plot(t_arr,1e3*waveforms[0],linewidth=linewidth,color='b')
                ax3 = fig.add_subplot(2,2,4)
                ax3.plot(t_arr,1e3*waveforms[1],linewidth=linewidth,color='r')
                for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
                    label.set_fontsize(tickSize)
                for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
                    label.set_fontsize(tickSize)
                ax2.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
                ax2.set_ylabel('Amplitude (mV)',fontsize=fontSize)
                ax2.set_title('$\sigma_x$ waveform',fontsize=fontSize)
                ax3.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
                ax3.set_title('$\sigma_z$ waveform',fontsize=fontSize)

                if tauk != 0:
                    textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$\n$\\nu$ = %.1f kHz\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],mu*1e3,AC_freq*1e-9,sigma*1e3,B0*1e3,tauk,nu,nAverages)
                elif nu == 0:
                    textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],mu*1e3,AC_freq*1e-9,sigma*1e3,B0*1e3,nu,nAverages)
                if plot_mode == 0:
                    # for plotting single realization
                    # background + generalized markovian data fit
                    if fitFunc == 'envelope':
                        ax1.plot(x_vector[1],decay(x_vector[1], fitted_pars[0], fitted_pars[1], fitted_pars[2]),'g',linewidth=linewidth,label='$\omega = 2\pi\\times%.1f$ MHz, $\\tau/\\tau_0$ = %.1f'%(B0*1e-3,fitted_pars[1][1]/fitted_pars[0][1]))
                    else:
                        ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
                    plt.gcf().text(0.925, 0.25, textstr, fontsize=fontSize,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
                    ax1.set_title('Ramsey Measurement %03d' %(iteration),size=fontSize)
                elif plot_mode == 1:
                    textstr = '$T_2^*$=%.2f $\mu$s\n$\sigma$ = %.1f mV'%(fitted_pars[3],sigma*1e3)
                    plt.gcf().text(0.925, 0.45, textstr, fontsize=fontSize+10,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
                elif plot_mode == 2:
                    textstr = '$\omega$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$'%(fitted_pars[1],fitted_pars[3],sigma*1e3,B0*1e3,nu)
                    plt.gcf().text(0.925, 0.45, textstr, fontsize=fontSize+10,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))

        elif (mu == 0 and B0 == 0 and sigma == 0) and simulation == 0: # plot data without any noise instances
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
        if fit_single_par_point == 1: # for plotting data averaged over many noise realizations post par sweep
            fig = plt.figure(figsize=((14,8)))
            ax1 = fig.add_subplot(111)
            # Plot data
            fontSize = 22
            tickSize = 20
            markersize= 10
            linewidth = 6

            # plot data
            ax1.plot(x_vector[0], y_vector[0], 'o', markersize = markersize, c='C0',label='$B_0 = 0$')
            if nu != 0:
                # if tauk/fitted_pars[0][1] > 100:
                label='$B_0$ = %.1f mV | $\\nu$ = %.1f kHz | $\\tau_k/\\tau_0$ = $\infty$'%(B0*1e3,nu)
                # else:
                    # label='$B_0$ = %.1f mV | $\\nu$ = %.1f kHz | $\\tau_k/\\tau_0$ = %.1f'%(B0*1e3,nu,tauk/fitted_pars[0][1])
            elif nu == 0:
                label='$B_0 = %.1f mV$ | $\\tau_k/\\tau_0 = %.3f$'%(B0*1e3,tauk/fitted_pars[0][1])
            ax1.plot(x_vector[1], y_vector[1], 'o', markersize = markersize, c='k',label=label)

            # plot fits
            ax1.plot(x_vector[0],decay(x_vector[0], fitted_pars[0][0], fitted_pars[0][1], fitted_pars[0][2]),'r',linewidth=linewidth,label='$\\tau_0 = %.1f \mu s$'%(fitted_pars[0][1]))
            ax1.plot(x_vector[1],decay(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], fitted_pars[1][2]),'g',linewidth=linewidth,label='$\\tau/\\tau_0$ = %.1f'%(fitted_pars[1][1]/fitted_pars[0][1]))

            ax1.legend(loc='upper right',prop={'size':22})
            ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
            ax1.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
            for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                label.set_fontsize(tickSize)
            if x_vector[0][-1] < x_vector[1][-1]:
                ax1.set_xlim([0,x_vector[0][-1]])
            else:
                ax1.set_xlim([0,x_vector[1][-1]])

        elif fit_single_par_point == 0 and simulation == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
            ax.set_ylabel('Digitizer Voltage (mV)')
            ax.set_xlabel('Pulse Separation ($\mu$s)')
            ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
            textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_2$=%.2f $\mu$s\n$\mu$ = %.3f V\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.3f V\n$B_0$ = %.3f V\n$\\tau_k$ = $\infty$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],mu,AC_freq*1e-9,sigma,B0,nAverages)
            ax.set_title('Echo Measurement %03d' %(iteration))
            plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

        elif simulation == 1:

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_vector, y_vector, 'bo', markersize = 3, label='Data')
            ax.set_ylabel('Fidelity')
            ax.set_xlabel('Time ($\mu$s)')
            ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r',label='Fit')
            ax.legend()
            ax.set_title('$B_0$ = %d kHz | $\\nu$ = %d kHz | $\\tau_k$ = %1.f $\mu$s'%(B0,nu,tauk))

    elif sequence == "T1":
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Delay ($\mu$s)')
        ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
        textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_1$=%.2f $\mu$s\n$\mu$ = %.3f V\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.3f V\n$B_0$ = %.2f V\n$\\tau_k$ = %.2f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],mu,AC_freq*1e-9,sigma,B0,nu,nAverages)
        ax.set_title('T1 Measurement %03d' %(iteration))
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)
    return fig

def plot_slice(B0,nu,tauk,sweep_var,ydata,error_T2,data='tau'):
    '''
    Plots a slice of a 3d plot

    Parameters
    ----------
    par1 : double
        B0
    par2 : double
        nu
    par3: double
        tau_k
    sweep_var: string
        which parameter goes on the xaxis
    ydata : double (1D array)
    data_type : string, optional
       Defines whether data are omega or tau. The default is 'tau'.
    '''

    fontsize = 14
    tickSize = 12
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt_options  = {
        'ecolor':      'r',
        'elinewidth':   5,
        'fmt':          '-o',
        'color':        'b',
        'capsize':      5,
        }

    if sweep_var == 'B0':
        ax.errorbar(B0,ydata,yerr=error_T2,**plt_options,label='$\\tau_k$ = %.1f $\mu$s'%(tauk[0]))
        ax.set_xlabel('$B_0$ (kHz)',fontsize=fontsize)
        ax.set_ylabel('$\\tau (\mu$s)',fontsize=fontsize)
    elif sweep_var == 'nu':
        ax.errorbar(nu,ydata,yerr=error_T2,**plt_options,label='$B_0$ = %d kHz\n$\\tau_k$ = %.1f $\mu s$'%(B0,tauk[0]))
        ax.set_xlabel('$\\nu$ (kHz)',fontsize=fontsize)
        ax.set_ylabel('$\\tau (\mu$s)',fontsize=fontsize)
    elif sweep_var == 'tau_k':
        ax.errorbar(tauk,ydata,yerr=error_T2,**plt_options,label='$B_0$ = %d kHz'%(B0))
        ax.set_xlabel('$\\tau_k$ ($\mu$s)',fontsize=fontsize)
        ax.set_ylabel('$\\tau (\mu$s)',fontsize=fontsize)
    else:
        raise Exception('Incorrect parameter combination!')

    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.legend(loc='center right',fontsize=fontsize-3)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(tickSize)
    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    plt.show()

def fit_sweep(par1,par2,par3,sweep='sweep_001',noiseType=1,nBackMeasurements=128,nMeasurements=128,sweep_type='nu_tauk',
              fileformat='new',fitFunc='',sequence='ramsey'):

    # initialize arrays
    T2_b_arr = np.zeros((len(par1),len(par2)))
    detun_b_arr = np.zeros((len(par1),len(par2)))
    T2_arr = np.zeros((len(par1),len(par2)))
    detun_arr = np.zeros((len(par1),len(par2)))
    error_T2_arr = np.zeros((len(par1),len(par2)))
    error_T2_b_arr = np.zeros((len(par1),len(par2)))

    if fitFunc == 'modDecayPlusCos':
        T2_2_arr = np.zeros((len(par1),len(par2)))
    else:
        T2_2_arr = []


    plt.close('all')
    for i in range(len(par2)):
        for j in range(len(par1)):
            if noiseType == 1:
                print('Now fitting B0 = %.1f mV | tau_k = %.3f us' %(par1[j]*1e3,par2[i]))
                fig,fitted_pars,error_T2,error_T2_b = plot_single_par_point(par1=par1[j], par2=par3[0], par3=par2[i],nBackMeasurements=nBackMeasurements,nMeasurements=nMeasurements, sweep=sweep,plot=1,noiseType=1,fileformat=fileformat,sequence=sequence,fitFunc=fitFunc)
                fig = plt.savefig(os.path.join('E:\\generalized-markovian-noise\\CandleQubit_6\\sweep_data\\%s\\'%sequence+sweep+'\\plot_images\\plot_B0_%d_uV_tau_%d_ns.png' %(round(par1[j]*1e6),round(par2[i]*1e3))) , bbox_inches='tight')

            elif noiseType == 2 and sweep_type == 'nu_tauk':
                # B0 = convertB0_kHz_to_mV(np.sqrt(2/9*1e6*(1-1/par2[i])**2+2*(2*np.pi*par1[j])**2)/(2*np.pi))*1e-3
                B0 = par3[0]
                print('Now fitting B0 = %.1f mV| nu = %.3f kHz| tau_k = %.3f us' %(B0*1e3,par1[j],par2[i]))
                fig,fitted_pars,error_T2,error_T2_b = plot_single_par_point(par1=B0, par2=par1[j], par3=par2[i],sweep=sweep,nBackMeasurements=nBackMeasurements,nMeasurements=nMeasurements, plot=1,noiseType=2,fileformat=fileformat,fitFunc=fitFunc,sequence=sequence)
                fig = plt.savefig(os.path.join('E:\\generalized-markovian-noise\\CandleQubit_6\\sweep_data\\%s\\'%sequence+sweep+'\\plot_images\\plot_B0_%d_uV_nu_%d_Hz_tau_%d_ns.png' %((par3[0]*1e6),round(par1[j]),round(par2[i]*1e3))) , bbox_inches='tight')

            elif noiseType == 2 and sweep_type == 'B0_nu':
                print('Now fitting B0 = %.1f mV| nu = %.3f kHz| tau_k = %.3f us' %(par1[j]*1e3,par2[i],par3[0]))
                fig,fitted_pars,error_T2,error_T2_b = plot_single_par_point(par1[j],par2[i],par3[0],sweep=sweep,nBackMeasurements=nBackMeasurements,nMeasurements=nMeasurements, plot=1,noiseType=2,fileformat=fileformat,fitFunc=fitFunc,sequence=sequence)
                fig = plt.savefig(os.path.join('E:\\generalized-markovian-noise\\CandleQubit_6\\sweep_data\\%s\\'%sequence+sweep+'\\plot_images\\plot_B0_%d_uV_nu_%d_Hz_tau_%d_ns.png' %(round(par1[j]*1e6),round(par2[i]),round(par3[0]*1e3))) , bbox_inches='tight')

            elif noiseType == 2 and sweep_type == 'B0_tauk':
                print('Now fitting B0 = %.1f mV| nu = %.3f kHz| tau_k = %.3f us' %(par1[j]*1e3,par3[0],par2[i]))
                fig,fitted_pars,error_T2,error_T2_b = plot_single_par_point(par1[j],par3[0],par2[i],sweep=sweep,nBackMeasurements=nBackMeasurements,nMeasurements=nMeasurements, plot=1,noiseType=2,fileformat=fileformat,fitFunc=fitFunc,sequence=sequence)
                fig = plt.savefig(os.path.join('E:\\generalized-markovian-noise\\CandleQubit_6\\sweep_data\\%s\\'%sequence+sweep+'\\plot_images\\plot_B0_%d_uV_nu_%d_Hz_tau_%d_ns.png' %(round(par1[j]*1e6),round(par3[0]),round(par2[i]*1e3))) , bbox_inches='tight')


            if noiseType == 1 and fitFunc != 'envelope':
                T2_arr[j,i] = fitted_pars[1][3]
                detun_arr[j,i] = fitted_pars[1][1]
            elif noiseType == 1 and fitFunc == 'envelope':
                T2_arr[j,i] = fitted_pars[1][1]
                detun_arr[j,i] = 0
                error_T2_arr[j,i] = error_T2
                T2_2_arr = []
            elif fitFunc == 'envelope' or sequence == 'echo':
                T2_arr[j,i] = fitted_pars[1][1]
                detun_arr[j,i] = 0
                error_T2_arr[j,i] = error_T2
                T2_2_arr = []
            elif fitFunc == 'modDecayPlusCos':
                T2_arr[j,i] = fitted_pars[1][3]
                T2_2_arr[j,i] = fitted_pars[1][5]
                error_T2_arr[j,i] = error_T2[1]
            elif fitFunc == 'modDampedCos':
                T2_arr[j,i] = fitted_pars[1][4]
                detun_arr[j,i] = fitted_pars[1][1]
                error_T2_arr[j,i] = error_T2

            T2_b_arr[j,i] = fitted_pars[0][1]
            error_T2_b_arr[j,i] = error_T2_b

    if fitFunc == 'modDecayPlusCos':
        T2_arr = [T2_arr,T2_2_arr]


    return detun_arr, T2_arr, T2_b_arr, error_T2_arr, error_T2_b_arr

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

def plot_single_par_point(par1,par2,par3,sweep,nMeasurements=100,nBackMeasurements=100,meas_device='CandleQubit_6',plot=1,noiseType=1,
                          fileformat='new',fitFunc='',sequence='ramsey',iteration=1):
    '''
    DESCRIPTION: Plots the averaged ramsey trace over the different noise realizations
    '''
    tdata_background,ydata_background, tdata, ydata, exp_pars = extract_data(sweep,par1, par2,par3,meas_device=meas_device,fileformat=fileformat,nBackMeasurements=nBackMeasurements,nMeasurements=nMeasurements, sequence=sequence)
    # tdata_background,ydata_background, tdata, ydata, exp_pars = extract_data(sweep,par1, par2,par3,meas_device=meas_device,fileformat=fileformat,nBackMeasurements=nBackMeasurements,nMeasurements=nMeasurements, sequence=sequence,iteration=iteration)

    #remove mean
    # for i in range(nMeasurements):
    #     ydata_background[i,:] = ydata_background[i,:] - np.mean(ydata_background[i,:])
    #     ydata[i,:] = ydata[i,:] - np.mean(ydata[i,:])

    # average
    ydata_avg_background = np.mean(ydata_background,axis=0)
    ydata_avg = np.mean(ydata,axis=0)

    #fit
    fitted_pars_background,error_b = fit_data(x_vector=tdata_background, y_vector=ydata_avg_background,sequence='echo',dt=tdata_background[-1]/len(tdata_background),verbose=0,**exp_pars,noiseType=1)
    fitted_pars,error = fit_data(x_vector=tdata, y_vector=ydata_avg,sequence=sequence,dt=tdata[-1]/len(tdata),**exp_pars,noiseType=noiseType,verbose=0,fitFunc=fitFunc)

    if fitFunc == 'modDampedCos':
        error_T2 = error[3]
    elif fitFunc == 'modDecayPlusCos':
        error_T2 = [error[3]]
        error_T2.append(error[6])
    elif fitFunc == 'envelope' or sequence == 'echo':
        error_T2 = error[1]
    elif noiseType == 1:
        error_T2 = error[3]

    #plot
    if plot == 1:
        fig = plot_data(awg=None,x_vector=[tdata_background,tdata],y_vector=[ydata_avg_background,ydata_avg],error=[error_b[1],error_T2],sequence=sequence,fitted_pars=[fitted_pars_background,fitted_pars],**exp_pars,fit_single_par_point=1,noiseType=noiseType,fitFunc=fitFunc)

    return fig,[fitted_pars_background,fitted_pars],error_T2,error_b[1]

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

def plot_single_shot(data_OFF,data_pi):
    states = []
    [states.append('|g>') for i in range(len(data_OFF.real))]
    [states.append('|e>') for i in range(len(data_pi.real))]
    data = {
            'I (mV)':   np.hstack((data_OFF.real*1e3,data_pi.real*1e3)),
            'Q (mV)':   np.hstack((data_OFF.imag*1e3,data_pi.imag*1e3)),
            'states':   states
                }
    dataF = pd.DataFrame(data=data)
    plot = sns.jointplot(data=dataF, x='I (mV)',y='Q (mV)',hue='states')
    # plot.fig.set_figwidth(7)
    # plot.fig.set_figheight(4)

def plot_noise(length=1000,gif_make=0,wfm1=[],wfm2=[]):
    mu = 325
    A_d = 200
    B0 = 100
    fontsize = 40
    linewidth = 10
    ticksize = 30
    if gif_make == 0:
        AC_stark_waveform = np.concatenate((0*np.ones(200),mu*np.ones(500),mu*np.ones(100),mu*np.ones(length)+np.random.normal(loc=0.325, scale=30, size=length),mu*np.ones(100),mu*np.ones(250),0*np.ones(200)))
        qubit_channel_waveform = np.concatenate((0*np.ones(200),A_d*np.zeros(500),A_d*np.ones(100),B0*gen_tel_noise(length, tau=0.75e6, dt=0.01),A_d*np.ones(100),0*np.ones(250),0*np.ones(200)))
    elif gif_make == 1:
        AC_stark_waveform = np.concatenate((wfm1[:800],wfm1[800:length+800],wfm1[-550:]))
        print(len(AC_stark_waveform))
        qubit_channel_waveform = np.concatenate((wfm2[:800],wfm2[800:length+800],wfm2[-550:]))
    t = np.linspace(0,1e-6,len(qubit_channel_waveform))
    fig = plt.figure(figsize=(20,16))
    ax = fig.add_subplot(111)
    ax.plot(t,AC_stark_waveform,c='r',linewidth=linewidth,label='AC Stark Channel Waveform ($\sigma_z$)')
    ax.set_xlabel('time ($\mu$s)',fontsize=fontsize)
    ax.plot(t,qubit_channel_waveform,linewidth=linewidth,label='Qubit Channel Waveform ($\sigma_x$)')
    ax.set_ylabel('Voltage (mV)',fontsize=fontsize)
    ax.set_ylim([-250,600])
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(ticksize)
    ax.legend(loc='upper center',prop={'size':30})
    plt.show()
    return fig,AC_stark_waveform,qubit_channel_waveform

def create_wfm_gif(Lmax=1000):
    numImag = 50
    filenames = []
    fig,AC_stark_waveform,qubit_channel_waveform = plot_noise(length=Lmax,gif_make=0)
    path = 'G:\Shared drives\LFL\Projects\Generalized Markovian noise\MarchMeeting2022\gif\\'
    for i in range(numImag):
        fig,wfm1,wfm2 = plot_noise(int(i*Lmax/numImag),gif_make=1,wfm1=AC_stark_waveform,wfm2=qubit_channel_waveform)
        filename = f'{i}.png'
        filenames.append(filename)
        plt.savefig('G:\Shared drives\LFL\Projects\Generalized Markovian noise\MarchMeeting2022\gif\\'+filename,  bbox_inches='tight')
        plt.close(fig)

   # build gif
    with imageio.get_writer(os.path.join(path,'ramsey_sequence.gif'), mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(os.path.join(path,filename))
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)

def gen_tel_noise(numPoints,tau,dt):

    signal = np.ones(numPoints)*(-1)**np.random.randint(0,2)
    for i in range(1,numPoints-1):
        if np.random.rand() < 1/(2*tau*1e-6/dt)*np.exp(-1/(2*tau*1e-6/dt)):
            signal[i+1] = - signal[i]
        else:
            signal[i+1] = signal[i]
    return signal


# def fit_single_instance(par1,par2,sweep='sweep_001',plot=1):
#     # start = time.time()
#     t_b, data_b , t , data = extract_data(sweep,par1=par1, par2 = par2)
#     # data = data[instance,:]
#     # detun , T2, error = pulse_plot1d(sequence='ramsey', dt=20/data.shape[1], plot=plot,x_vector=t, y_vector=np.mean(data,axis=0),AC_pars=[0.6,0.08],RT_pars=[par1,par2])
#     fit_beats(sequence='ramsey', dt=5/data.shape[1], plot=plot,x_vector=t, y_vector=-np.mean(data,axis=0),AC_pars=[0.6,0.08],RT_pars=[par1,par2])
#     end = time.time()

#     # return detun, T2

# def extract_data(sweep,par1,par2,par3,meas_device='CandleQubit_6'):
#     start = time.time()
#     filename = 'B0_%d_uV_nu_%d_kHz_tau_%d_ns' %(round(par1*1e6),round(par3*1e3),round(par2*1e3))
#     # filename = 'RTN_B0_%d_uV_tau_%d_ns' %(round(par1*1e6),round(par2*1e3))
#     with open("E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep,filename)) as datafile:
#         csv_reader = csv.reader(datafile,delimiter=',')
#         file_data = list(csv_reader)

#         for row in file_data:
#             if row[0] == "Background Time Data":
#                 print(file_data[file_data.index(row)+1])
#                 tdata_background = np.array(file_data[file_data.index(row)+1],dtype=float)
#             if row[0] == "Time Data":
#                 datafile.seek(0)
#                 tdata= np.array(file_data[file_data.index(row)+1],dtype=float)
#             if row[0] == "Background Data: Channel 1":
#                 line_start_background = file_data.index(row) + 1
#             if row[0] == "Background Data: Channel 2":
#                 line_end_background = file_data.index(row) - 1
#             if row[0] == "Data: Channel 1":
#                 line_start = file_data.index(row) + 1
#             if row[0] == "Data: Channel 2":
#                 line_end = file_data.index(row) - 1


#         datafile.seek(0)
#         # extract traces
#         ydata_background = np.zeros((line_end_background-line_start_background+1,len(tdata_background)))
#         line = line_start_background
#         while line <= line_end_background:
#             datafile.seek(0)
#             # print(np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32))
#             ydata_background[line-line_start_background,:] = np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32)
#             line += 1
#             # print(ydata_background)

#         datafile.seek(0)
#         trace = np.array(next(itertools.islice(csv_reader, line_start,None)),dtype=np.float32)
#         # ydata = np.zeros((line_end-line_start+1,len(trace)))
#         ydata = np.zeros((line_end-line_start+1,len(tdata)))
#         line = line_start
#         while line <= line_end:
#             datafile.seek(0)
#             # data = np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32)
#             # ydata[line-line_start,:] = data[:148]
#             ydata[line-line_start,:] = np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32)[:]
#             line += 1

#         #extract data point parameters (qubitDriveFreq, ACstarkFreq, etc.)
#         datafile.seek(0)
#         keys = file_data[0]
#         values =  file_data[1]
#         dictionary = dict.fromkeys(keys,0)

#         i = 0
#         for key in dictionary:
#             if key == 'AC_pars' or key == 'RT_pars':
#                 pars = np.zeros((3,1))
#                 values[i] = values[i].replace('[','')
#                 values[i] = values[i].replace(']','')
#                 values[i] = values[i].replace(',','')
#                 pars[0] = float(values[i][1:values[i].find(' ')])
#                 pars[1] = float(values[i][values[i].find(' ')+1:values[i].find(' ',values[i].find(' ')+1)])
#                 pars[2] = float(values[i][values[i].find(' ',values[i].find(' ')+1):])
#                 dictionary[key] = pars
#             try:
#                 dictionary[key] = float(values[i])
#             except:
#                 dictionary[key] = values[i]
#             i += 1

#     end = time.time()
#     print(dictionary)
#     print('Took %d seconds to extract data'%(end-start))
#     return tdata_background, ydata_background, tdata, ydata, dictionary
# def plot_single_par_point(par1,par2,par3,sweep,meas_device='CandleQubit_6',plot=1):
#     '''
#     DESCRIPTION: Plots the averaged ramsey trace over the different noise realizations
#     '''
#     tdata_background,ydata_background, tdata, ydata, exp_pars = extract_data(sweep, par1, par2,par3,meas_device=meas_device)
#     # average
#     ydata_avg_background = np.mean(ydata_background,axis=0)
#     ydata_avg = np.mean(ydata,axis=0)
#     #fit
#     fitted_pars,error = fit_data(x_vector=tdata, y_vector=ydata_avg,dt=tdata[-1]/len(tdata),**exp_pars)
#     fitted_pars_background,error_b = fit_data(x_vector=tdata_background, y_vector=ydata_avg_background,dt=tdata_background[-1]/len(tdata_background),**exp_pars)
#     #plot
#     if plot == 1:
#         fig = plot_data(x_vector=[tdata_background,tdata],y_vector=[ydata_avg_background,ydata_avg],fitted_pars=[fitted_pars_background,fitted_pars],**exp_pars,fit_single_par_point=1)

#     return fig,[fitted_pars_background,fitted_pars]

# def fit_beats(sequence,x_vector,y_vector,dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,AC_pars=[0,0],RT_pars=[0,0],pi2Width=50e-9,fitting = 1, plot=1,save_fig = 0,iteration=1):
#     x_vector = x_vector*1e6
#     abs_camplitude = np.abs(y_vector*1e3)
#     amp = abs_camplitude[0]-abs_camplitude[-1]
#     offset = np.mean(abs_camplitude)
#     f1 = 15
#     f2 = 1
#     tau = 15
#     phi1 = 0
#     phi2 = 0
#     p0 = [amp,f1,f2,phi1,phi2,tau,offset]
#     lb = [-1000,-10,-10,-np.pi,-np.pi,0,-np.inf]
#     ub = [1000,20,20,np.pi,np.pi,30,np.inf]
#     fitted_pars, covar = scy.optimize.curve_fit(beats, x_vector, abs_camplitude,p0=p0,bounds=(lb,ub),xtol=1e-12,maxfev=6000)
#     f1 = best_vals[1]
#     f2 = best_vals[2]
#     T_phi = best_vals[5]
#     error = np.sqrt(abs(np.diag(covar)))
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.plot(x_vector, abs_camplitude, '-o', markersize = 3, c='C0')
#     ax.set_ylabel('Digitizer Voltage (mV)')
#     ax.set_xlabel('Pulse Separation ($\mu$s)')
#     ax.plot(x_vector,beats(x_vector,best_vals[0],best_vals[1],best_vals[2],best_vals[3],best_vals[4],best_vals[5],best_vals[6]),'r')
#     textstr = '$A$ = %.2f mV\n$\omega_1$=%.2f MHz\n$\omega_2$=%.2f MHz\n$T_2^*$=%.2f $\mu$s\n$B_0$ = %.2f V\n$\\tau_k$ = %.2f $\mu s$'%(best_vals[0],best_vals[1],best_vals[2],best_vals[5],RT_pars[0],RT_pars[1])
#     ax.set_title('Ramsey Measurement %03d' %(iteration))
#     plt.gcf().text(1, 0.25, textstr, fontsize=14)
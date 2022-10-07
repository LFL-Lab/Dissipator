# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:46:39 2022

@author: haimeng zhang <haimeng@usc.edu>
"""

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Basic functions for fitting
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def expCosine(x_data, a, freq, t2, phase, offset):
    return (a*np.exp(-x_data/t2)*np.cos(2*np.pi*(freq)*x_data + phase) + offset)
keys1 = ['a', 'freq', 't2', 'phase', 'offset']

def exp(x_data, a, t2, offset):
     return (a*np.exp(-x_data/t2)+offset)
keys2 = ['a', 't2', 'offset']

def expCosCos(x_data, a, freq1, freq2, t2, phase1, phase2, offset):
    return (a*np.exp(-x_data/t2)*np.cos(2*np.pi*(freq1)*x_data + phase1)*np.cos(2*np.pi*(freq2)*x_data + phase2) + offset)
keys3 = ['a', 'freq1', 'freq2', 't2', 'phase1', 'phase2', 'offset']

def fft_freq(ydata,xdata,plot=True, **kwargs):
    """
    FFT of ydata, return signal period in ns
    d: sample distance, 1/sampling rate
    n: optional, length of data by defalut
    test with: fft_period(ydata,xdata,n=300)
    """
    if 'height' in kwargs:
        height = kwargs.get('height')
    else:
        height = 10
    n = len(xdata)
    d =(xdata[1]-xdata[0]) # unit: us
    fydata = np.abs(np.fft.fft(ydata,n))
    freq = np.fft.fftfreq(n,d)
    nf = len(freq)
    peaks, _ = find_peaks(fydata, height=height)
    guess_freq = freq[peaks]
    if plot:
        plt.plot(freq[:nf//2+1],fydata[:nf//2+1:])
        plt.plot(guess_freq, fydata[peaks], "x")
    for f in guess_freq:
        print("The guess freq {} MHz".format(f))
    return guess_freq


def perform_fit(xdata, ydata, function=expCosine, keys=keys1, phase_0=0, plot=False, **kwargs):
    if 't2_0' in kwargs:
        t2_0 = kwargs.get('t2_0')
    else:
        t2_0 = xdata[-1]
    amp_0 = abs(max(ydata) - min(ydata))
    if function.__name__ == 'expCosine':
        offset_0 = np.mean(ydata)
        if 'freq_0' in kwargs:
            freq_0 = kwargs.get('freq_0')
        else:
            freq_0 = 2 / xdata[-1]
        guess = [amp_0 / 2, freq_0, t2_0, phase_0, offset_0]

    if function.__name__ == 'expCosCos':
        offset_0 = np.mean(ydata)
        if 'freq1_0' in kwargs:
            freq1_0 = kwargs.get('freq1_0')
        else:
            freq1_0 = 2 / xdata[-1]
        if 'freq2_0' in kwargs:
            freq2_0 = kwargs.get('freq2_0')
        else:
            freq2_0 = 0
        if 'phase1_0' in kwargs:
            phase1_0 = kwargs.get('phase1_0')
        else:
            phase1_0 = 0
        if 'phase2_0' in kwargs:
            phase2_0 = kwargs.get('phase2_0')
        else:
            phase2_0 = 0
        guess = [amp_0, freq1_0, freq2_0, t2_0, phase1_0, phase2_0, offset_0]


    elif function.__name__ == 'exp':
        offset_0 = np.min(ydata)
        guess = [amp_0 / 2, t2_0, offset_0]

    pars, covar = curve_fit(function, xdata, ydata, p0=guess)
    error = np.sqrt(abs(np.diag(covar)))
    perr = np.sqrt(np.diag(covar))
    fit_result = dict(zip(keys, pars))
    fits_are_good = (abs(error / pars))[1] < 1

    if plot:
        plt.plot(xdata, ydata, 'x', label='data')

        if function.__name__ == 'expCosine':
            yFit = expCosine(xdata, fit_result['a'], fit_result['freq'],
                             fit_result['t2'], fit_result['phase'], fit_result['offset'])
            label = 'fit: ' + r'$%.1fe^{-t/%.1f}\cos(2\pi*%.1f*t+%.1f)+%.1f$' % (fit_result['a'], fit_result['t2'],
                                                                                 fit_result['freq'],
                                                                                 fit_result['phase'],
                                                                                 fit_result['offset'])

        elif function.__name__ == 'exp':
            yFit = exp(xdata, fit_result['a'], fit_result['t2'], fit_result['offset'])
            label = 'fit: ' + r'$%.1fe^{-t/%.1f}+%.1f$' % (fit_result['a'], fit_result['t2'], fit_result['offset'])

        elif function.__name__ == 'expCosCos':
            yFit = expCosCos(xdata, fit_result['a'], fit_result['freq1'], fit_result['freq2'],
                             fit_result['t2'], fit_result['phase1'], fit_result['phase2'], fit_result['offset'])
            label = 'fit: ' + r'$%.1fe^{-t/%.1f}\cos(2\pi*%.1f*t+%.1f)\cos(2\pi*%.1f*t+%.1f)+%.1f$' % (
            fit_result['a'], fit_result['t2'],
            fit_result['freq1'], fit_result['phase1'],
            fit_result['freq2'], fit_result['phase2'],
            fit_result['offset'])

        plt.plot(xdata, yFit, '-',
                 label=label)
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    return fit_result
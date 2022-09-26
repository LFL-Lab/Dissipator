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
from meas_utilities import *
from collections import OrderedDict

##### Functions for fitting #####

def expCosine(xdata, a, freq, t2, phase, offset):
    return (a * np.exp(-xdata / t2) * np.cos(2*np.pi * freq *xdata + phase) + offset)

def Exp(xdata, a, t2, offset):
     return (a * np.exp(-xdata / t2) + offset)

def expCosCos(xdata, a, freq1, freq2, t2, phase1, phase2, offset):
    return (a * np.exp(-xdata / t2) * np.cos(2*np.pi * freq1 * xdata + phase1) * np.cos(2*np.pi * freq2 * xdata + phase2) + offset)

"""Default values for parameters of each function. Since order matters when inputting
the guess to the fitting utility, these are stored as a Dict of OrderedDicts.
"""
def guess(xdata, ydata):
    
    return {Exp       : OrderedDict(
                        {'a'      : abs(max(ydata) - min(ydata)) / 2,
                         't2'     : xdata[-1],
                         'offset' : np.min(ydata)
                        }),
            
            expCosine : OrderedDict(
                        {'a'     : abs(max(ydata) - min(ydata)) / 2,
                         'freq'  : 2 / xdata[-1],
                         't2'    : xdata[-1],
                         'phase' : 0,
                         'offset': np.mean(ydata)
                        }),
                
            expCosCos : OrderedDict(
                        {'a'      : abs(max(ydata) - min(ydata)),
                         'freq1'  : 2 / xdata[-1],
                         'freq2'  : 0,
                         't2'     : xdata[-1],
                         'phase1' : 0,
                         'phase2' : 0,
                         'offset' : np.mean(ydata)
                        })
            }

"""Given kwargs input to perform_fit, modify the default OrderedDict returned by 'guess'
with the values specified in kwargs.
"""
def getguess(fitfunc, xdata, ydata, kwargs):
    
    # get default guess corresponding to fitfunc
    g = guess(xdata, ydata)[fitfunc]
    
    # extract keys / values in kwargs that correspond to default values for fitting 
    subkwargs = subkeywords(g.keys(), kwargs)
    
    # replace with any input keyword argument values
    for (key, value) in subkwargs.items():
        g[key] = value
    
    return g
        
"""
Return a new OrderedDict that is the subset of items in 'kwargs' that have keys in 'keys';
these items are removed from kwargs.
"""
def subkeywords(keys, kwargs):
    
    supplied_kwargs = OrderedDict()
    for key in keys:
        if key in kwargs:
            supplied_kwargs[key] = kwargs.pop(key) 
            
    return supplied_kwargs


"""
Fit 'xdata' and 'ydata' using function 'fitfunc'. 
Optionally, plot the data and the fit result.
Returns the fit result as a dict of parameter names and values; also returns the plot if applicable.
"""
def perform_fit(fitfunc, xdata, ydata,  plot = False,       # plot the data and fit result
                                        maketitle = True,   # give the plot a title that prints the fit result
                                        title = "",         # title string for plot (prepend to the fit result title, if applicable)
                                        precision = 4,      # scientific precision of fit parameters printed in title
                                        **kwargs):          # non-default values for initial guess parameters, as well as plotting keyword arguments

    # get guess for fitting, extract names of fit parameters as 'keys', and default values as 'p0'
    g = getguess(fitfunc, xdata, ydata, kwargs)
    keys, p0 = list(g.keys()), list(g.values())

    # fit
    pars, covar = curve_fit(fitfunc, xdata, ydata, p0=p0)
    error = np.sqrt(abs(np.diag(covar)))
    perr = np.sqrt(np.diag(covar))
    fit_result = dict(zip(keys, pars))
    fits_are_good = (abs(error / pars))[1] < 1
    
    
    if plot:
        # clear plotting palette
        plt.clf()
        
        # make customized title based on fit results
        if maketitle:

            title += '\n' + make_title(paramdict = fit_result, 
                                precision = precision, 
                                scientific = True)
            
        # use subplots so that plotting kwargs can be automatically processed using plt.setp
        fig, ax = plt.subplots()
        ax.plot(xdata, ydata, 'x', label='data')
        plt.setp(ax, **kwargs)
        plt.setp(ax, title = title)
        
        # evaluate the fitting function at different xdata values, based on fit results
        yfit = fitfunc(xdata, *fit_result.values())

        # add fit line and label to plot
        ax.plot(xdata, yfit, '-', label = "fit")
        ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        fig.tight_layout()
        
        return fit_result, fig
        
    return fit_result




def fft_freq(ydata, xdata, plot=True, **kwargs):
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



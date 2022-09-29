import matplotlib.pyplot as plt
import numpy as np
import os
import csv

import seaborn as sns
import pandas as pd

from config import *


# where figures are saved
datapath = f'D:\\Users\\lfl\\data\\candle\\res_{resFreq/1e9}GHz\\'
foldername = datapath + 'qubit_spec_figures\\'


def heatplot(xaxis, yaxis, data, xlabel = "", ylabel = "", normalize=False, colorbarlabel = 'log mag', **kwargs):
    
    if normalize:
        colorbarlabel += ' (normalized)'
        
    df = pd.DataFrame(data, columns = xaxis, index = yaxis)
    
    if normalize:
        df = df.apply(lambda x: (x-x.mean())/x.std(), axis = 1)
    
    hm = sns.heatmap(df, cmap = 'viridis', cbar_kws={'label': colorbarlabel}, **kwargs)
    hm.set_xlabel(xlabel, fontsize=12)
    hm.set_ylabel(ylabel, fontsize=12)
    
    plt.tight_layout()
     
    return hm;
    
#I, Q, freqs, attens = scan_resonator_spec(attens = np.arange(50, 0, -1));

class IQData(object):
    pass



"Not yet working..."
def plotIQs(*datas: IQData, title = "spectroscopy", minfreq = None, maxfreq = None):
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    fig.suptitle(title, fontsize=11)
    
    ax1.set_ylabel("mag (V)")
    ax2.set_xlabel("freq (GHz)")
    ax2.set_ylabel("phase (rad)")
       
    
    for data in datas:
        start, stop = startstop(minfreq, maxfreq, data.freqs)
        
        I = data.I[start:stop]
        Q = data.Q[start:stop]
        freqs = data.freqs[start:stop]
        
        ax1.plot(freqs, mag(I, Q), alpha=0.7)
        ax2.plot(freqs, mag(I, Q), alpha=0.7)
       
    return fig;
    


def plotIQ(arg1, *args, **kwargs):
    
    if isinstance(arg1, np.ndarray):
        return plotIQ_both(arg1, *args, **kwargs)
    else:
        return plotIQ_data(arg1, *args, **kwargs)
        

def startstop(minfreq, maxfreq, freqs):
    
    if minfreq is None:
        start = 0
    else:
        start = np.where(freqs == minfreq)[0][0]
        
    if maxfreq is None:
        stop = len(freqs) - 1
    else:
        stop = np.where(freqs == maxfreq)[0][0]
        
    return start, stop
    

"""
plotIQ_both(
        I,      :: list of averaged I quadrature values, as function of frequency
        Q,      :: list of averaged Q quadrature values, as function of frequency
        freqs,  :: absolute frequencies
        title = "spectroscopy"):

Returns: Plot of magnitude and phase of the signal as a function of frequency.
"""
def plotIQ_both(I, Q, freqs, title = "spectroscopy", minfreq = None, maxfreq = None):
    
    start, stop = startstop(minfreq, maxfreq, freqs)   
    
    I = I[start:stop]
    Q = Q[start:stop]
    freqs = freqs[start:stop]
    
    mag = np.sqrt(I**2 + Q**2);
    phase = np.unwrap(np.angle(I + 1j*Q))
        
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    fig.suptitle(title, fontsize=11)
    
    ax1.plot(freqs, mag)
    ax1.set_ylabel("mag (V)")
    
    ax2.plot(freqs, phase)
    ax2.set_xlabel("freq (GHz)")
    ax2.set_ylabel("phase (rad)")
    
    plt.tight_layout()
    
    return fig;

def plotIQ_new(xdata, I, Q, title = "", xmin = None, xmax = None, xlabel = ""):
    
    start, stop = startstop(xmin, xmax, xdata)   
    
    I = I[start:stop]
    Q = Q[start:stop]
    xdata = xdata[start:stop]
        
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    fig.suptitle(title, fontsize=11)
    
    ax1.plot(xdata, mag(I, Q))
    ax1.set_ylabel("mag (V)")
    
    ax2.plot(xdata, phase(I, Q))
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("phase (rad)")
    
    plt.tight_layout()
    
    return fig;

def plotIQ_quads(xdata, I, Q, title = "", xmin = None, xmax = None, xlabel = ""):
    
    start, stop = startstop(xmin, xmax, xdata)   
    
    I = I[start:stop]
    Q = Q[start:stop]
    xdata = xdata[start:stop]
        
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    fig.suptitle(title, fontsize=11)
    
    ax1.plot(xdata, I)
    ax1.set_ylabel("I (V)")
    
    ax2.plot(xdata, Q)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Q (V)")
    
    plt.tight_layout()
    
    return fig;



def plotIQ_data(data: IQData, title = "spectroscopy", minfreq = None, maxfreq = None):
    
    return plotIQ_both(data.I, data.Q, data.freqs, title=title, minfreq=minfreq, maxfreq=maxfreq);


mag = lambda I, Q : np.sqrt(I**2 + Q**2)
phase = lambda I, Q : np.unwrap(np.angle(I + 1j*Q))

def mags(allI, allQ):
    return np.array([mag(I, Q) for (I, Q) in zip(allI, allQ)])

def import_qubit_spec_scan(csvfilename):
    
    with open(csvfilename, 'r') as f:
        
        reader = csv.reader(f)
        freqs = np.array([float(i) for i in next(reader)])
        I = np.array([float(i) for i in next(reader)])
        Q = np.array([float(i) for i in next(reader)])
        
    return freqs, I, Q


def plot_sweep_data_overlaid(foldername = foldername + "FF_sweep\\", subfoldername = lambda val: f'DC_bias_{val}_uA',
                           sweepvals = [0.0, 100.0, 200.0, 300.0, 400.0]):
    
    maglist = []
    phaselist = []
    datas = []
    
    for val in sweepvals:
        
        data = IQData()
        
        subfolder = subfoldername(val)
        path = foldername + "/" + subfolder
        files = os.listdir(path)
        
        for filename in files:
            
            if '.csv' in filename:
                csvfile = filename
                break
        
        data.freqs, data.I, data.Q = import_qubit_spec_scan(path + '/' + csvfile)
        datas.append(data)
        
        
    return plotIQs(*datas)

 
            
"""
Plots sweep data generated by sweep_FF. Takes as argument foldername, the same as input to sweep_FF.

COMMENT: This was meant to be a heatplot but it needs to be modified. Also, note that the plot will only make sense
visually if there are a similar number of datapoints in the flux dimension as in the frequency dimension.
"""            
def plot_sweep_data(foldername, subfoldername = lambda val: f'DC_bias_{val}_uA',
                           sweepvals = [0.0, 100.0, 200.0, 300.0, 400.0]):
    
    maglist = []
    phaselist = []
    
    for val in sweepvals:
        
        subfolder = subfoldername(val)
        path = foldername + "/" + subfolder
        files = os.listdir(path)
        
        for filename in files:
            
            if '.csv' in filename:
                csvfile = filename
                break
        
        freqs, I, Q = import_qubit_spec_scan(path + '/' + csvfile)
        maglist.append(mag(I, Q))
        phaselist.append(phase(I, Q))
        
        
    # make 3D plot
    x_var = freqs
    y_var = sweepvals    
    X, Y = np.meshgrid(x_var, y_var)
    
    
    fig, ax = plt.subplots(1,1)
    cp = ax.contourf(X * 1e-9, Y, maglist)
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('qubit spectroscopy')
    ax.set_xlabel('frequency (GHz)')
    ax.set_ylabel('bias current (uA)')
    plt.show()
        
    
    
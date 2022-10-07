import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import scipy.optimize
import pandas as pd
import os

from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager

from slacker import sendslack
from datetime import datetime, timedelta


def gettimedelta(progress_bar):
    
    elapsed = progress_bar.format_dict["elapsed"]
    rate = progress_bar.format_dict["rate"]
    seconds_remaining = (progress_bar.total - progress_bar.n) / rate
    delta = timedelta(seconds = seconds_remaining)
    now = datetime.now()
    later = now + delta
    return later.strftime("%H:%M:%S")

def datestring():
    
    now = datetime.now()
    string = now.strftime("%D")
    return string.replace("/", "-")


def timestring():
    
    now = datetime.now()
    string = now.strftime("%D_%H%M%S")
    return string.replace("/", "-")



def make_name(  s1 = "_",
                s2 = "_",
                paramdict = {}, 
                precision = 2, 
                scientific = False, 
                **kwargs):
    
   name = ""
   
   for param, value in paramdict.items():
       if scientific:
           value = np.format_float_scientific(value, precision = precision, trim = 'k')
       name += param + s1 + str(value) + s2
   
   for param, value in kwargs.items():
       if scientific:
           value = np.format_float_scientific(value, precision = precision, trim = 'k')
       name += param + s1 + str(value) + s2
       
   return name[:-1]


def make_title(*args, **kwargs):   return make_name(s1 = " = ", s2 = ", ", *args, **kwargs)[:-1]  

def writedata(name, datadict = {}, **kwargs):
    
    if len(datadict) > 0 and len(kwargs.items()) > 0:
        raise Exception("Please input EITHER a dictionary with keys as headers and values as data arrays, OR use keyword arguments with symbols as headers and values as data arrays.")
        
    if len(kwargs.items()) > 0:
        datadict = kwargs.items()
      
    df = pd.DataFrame(datadict)
    df.to_csv(name, index = False)
    
    
def getdata(name): return pd.read_csv(name)
    
    
def makefolder(datapath, name):
    
    foldername = datapath + name + "_" + datestring() + "\\"
    if not os.path.isdir(foldername):
        os.mkdir(foldername);
    
    return foldername

    

# https://stackoverflow.com/questions/16716302/how-do-i-fit-a-sine-curve-to-my-data-with-pylab-and-numpy
def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = np.array(tt)
    yy = np.array(yy)
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}


def fit_exp(tt, yy):

    def expfunc(x, a, b, c):  return a * np.exp(-b * x) + c
    
    guess_amp = np.std(yy) * 2.**0.5
    guess_rate = 1/ tt[-1]   # excluding the zero frequency "peak", which is related to offset
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, guess_rate, guess_offset])
    
    popt, pcov = scipy.optimize.curve_fit(expfunc, tt, yy, p0=guess)
    a, b, c = popt
    fitfunc = lambda t: a * np.exp(-b * t) + c
    return {"amp": a, "offset": c, "gamma": b, "T1": 1/(2*b), "fitfunc": fitfunc, "maxcov": np.max(pcov)}
    




"""
getIQ(
      config,       :: configuration dictionary (usually imported from a file named config.py),
                       which is interpretted by QUA to determine settings on OPX 
      jobtype):     :: QUA jobtype (e.g. qubit_spec, rr_spec)
"""
def getIQ(config, jobtype, showprogress = False, n_avg = 1000, notify = False):
    
    
    # Open Communication with the Server #
    qmm = QuantumMachinesManager()

    qm = qmm.open_qm(config)
    job = qm.execute(jobtype)
    res_handles = job.result_handles
     
    # res_handles.wait_for_all_values()

    I_handle = res_handles.get("I")
    Q_handle = res_handles.get("Q")
    I_handle.wait_for_values(1)
    Q_handle.wait_for_values(1)
    
    if showprogress:
        n_handle = res_handles.get("n")
        n0 = 0
        with tqdm(total = n_avg) as progress_bar:
            while(I_handle.is_processing()):
                
                n = n_handle.fetch_all()
                
                if notify and n0 <= 100 and n > 100:
                    
                    later = gettimedelta(progress_bar)
                    sendslack(message = f":zap: Measurement running. The measurement will be complete at approximately " + later + ".")
                
                Δn = n - n0
    
                if Δn > 0:
                    progress_bar.update(Δn)
                    n0 = n
                    
               
                    
          
    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    
    qmm.close_all_quantum_machines() # 
    
    return I, Q;
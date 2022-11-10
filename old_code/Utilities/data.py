import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import scipy.optimize
import pandas as pd
import os
import platform

from datetime import datetime, timedelta


# given a tqdm object progress_bar, gives an estimate in hours, minutes, and seconds
# of the total time remaining in the process
def gettimedelta(progress_bar):
    
    elapsed = progress_bar.format_dict["elapsed"]
    rate = progress_bar.format_dict["rate"]
    seconds_remaining = (progress_bar.total - progress_bar.n) / rate
    delta = timedelta(seconds = seconds_remaining)
    now = datetime.now()
    later = now + delta
    return later.strftime("%H:%M:%S")


# create a datestamp 
def datestring():
    
    now = datetime.now()
    string = now.strftime("%D")
    return string.replace("/", "-")


# create a timestamp from the date, hour minute, and second
def timestring():
    
    now = datetime.now()
    string = now.strftime("%D_%H%M%S")
    return string.replace("/", "-")


# given keys and values input as a dict 'paramdict' or as 'kwargs', generate
# a string encoding the keys and values in an acceptable filename
def make_name(  s1 = "_",           # separator between key and its value
                s2 = "_",           # separator between value and next key
                paramdict = {},     # keys and values to put in name
                precision = 2,      # scientific precision of values printed
                scientific = False, # use scientific notation
                **kwargs):          # keys and values to put in name
    
   # merge the two input dictionaries
   parameters = {**paramdict, **kwargs} 
   
   # loops over keys and values in parameters and generate string
   names = ""
   for name, value in parameters.items():
        if scientific:
            value = np.format_float_scientific(value, precision = precision, trim = 'k')
        names += name + s1 + str(value) + s2
    
   # return names, cutting off last '_' from string
   return names[:-1]


# creates a plot title similar to the filename returned by 'make_name'
def make_title(*args, **kwargs):   
    return make_name(s1 = " = ", s2 = ", ", *args, **kwargs)[:-1]  


def writedata(name, datadict = {}, **kwargs):
    
    # merge kwargs into datadict
    datadict.update(kwargs)
      
    # create dataframe and write to csv
    df = pd.DataFrame(datadict)
    df.to_csv(name, index = False)
    
# import data, in dataframe format, from csv file, given the file 'name' to import
def getdata(path, *paths, varnames = []): 
    
    # read in dataframe from csv
    data = pd.read_csv(getpath(path, *paths))
    
    # if user inputs variable names, return only those variables as a list of numpy arrays
    # otherwise, return the dataframe
    if varnames:
        return map(lambda v : np.array(data[v]), varnames)
        
    else:
        return data


# join elements of paths together and return in appropriate format for the os (Windows, Linux, or Mac)
# (on Windows, 'normpath' converts '/' to '\')
def getpath(path, *paths):
    newpath = os.path.join(path, *paths)
    return os.path.normcase(newpath)
    
def Volt2dBm(data):
    '''
    converts from voltage into dBm
    '''
    return 10*np.log10(1e3*data**2/50)

def Watt2dBm(x):
    '''
    converts from units of Watts to dBm
    '''
    return 10.*np.log10(x*1000.)

# make folder 'name' at location 'datapath' (cross-platform)
def makefolder(path, *paths, date = True):
    
     # create the normalized pathname 
    foldername = getpath(path, *paths)
    
    if date:
        foldername += "_" + datestring()
    
    # create the directory if it doesn't exist
    # if not os.path.isdir(foldername):
    #     os.mkdir(foldername);
    if not os.path.isdir(foldername):
        os.makedirs(foldername)
   
    return foldername



def clk(time_in_ns):
    
    if not time_in_ns % 4 == 0:
        raise Exception("Please enter a time in ns as a multiple of 4 for conversion to clock cycles.")
        
    if not type(time_in_ns) is int:
        raise Exception("Please enter a time in ns as an integer for conversion to clock cycles.")
    
    return time_in_ns // 4

    
# fitting functions moved to fitting.py, but still need to replace code in rabi

# # https://stackoverflow.com/questions/16716302/how-do-i-fit-a-sine-curve-to-my-data-with-pylab-and-numpy
# def fit_sin(tt, yy):
#     '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
#     tt = np.array(tt)
#     yy = np.array(yy)
#     ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
#     Fyy = abs(np.fft.fft(yy))
#     guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
#     guess_amp = np.std(yy) * 2.**0.5
#     guess_offset = np.mean(yy)
#     guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

#     def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
#     popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
#     A, w, p, c = popt
#     f = w/(2.*np.pi)
#     fitfunc = lambda t: A * np.sin(w*t + p) + c
#     return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}


# def fit_exp(tt, yy):

#     def expfunc(x, a, b, c):  return a * np.exp(-b * x) + c
    
#     guess_amp = np.std(yy) * 2.**0.5
#     guess_rate = 1/ tt[-1]   # excluding the zero frequency "peak", which is related to offset
#     guess_offset = np.mean(yy)
#     guess = np.array([guess_amp, guess_rate, guess_offset])
    
#     popt, pcov = scipy.optimize.curve_fit(expfunc, tt, yy, p0=guess)
#     a, b, c = popt
#     fitfunc = lambda t: a * np.exp(-b * t) + c
#     return {"amp": a, "offset": c, "gamma": b, "T1": 1/(2*b), "fitfunc": fitfunc, "maxcov": np.max(pcov)}
    



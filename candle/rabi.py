import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
import csv

from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager

from config import *
from plotting import *
from meas_utilities import *

import sys
sys.path.append("D:\Program Files\Keysight\Labber\Script")
import Labber

# where to save figures
datapath = f'D:\\Users\\lfl\\data\\candle\\res_{resFreq/1e9}GHz\\'
if not os.path.isdir(datapath):
    os.mkdir(datapath);
foldername = datapath + 'rabi_measurements' + datestring()
if not os.path.isdir(foldername):
    os.mkdir(foldername);



def timerabi(t0 = 40,         # minimum pulse duration in nanoseconds
              tf = 10e3,    # maximum pulse duration in nanoseconds
              dt = 500,        # step of sweep in nanoseconds
              n_avg = 1000,  # number of averages
              amp_q_scaling = 1,
              makeplots = True,
              magphase = False,
              foldername = foldername + "\\time_rabi"):
    
    if not os.path.isdir(foldername):
        os.mkdir(foldername);
        os.mkdir(foldername + "\\figures")
        os.mkdir(foldername + "\\data")

    t_min = round(t0 / 4)  # minimum pulse duration in clockcycles 
    t_max = round(tf / 4)  # maximum pulse duration in clockcycles
    step = round(dt / 4)   # step of sweep in clockcycles
              
    times = np.arange(t_min, t_max + step/2, step, dtype=int) # in clockcycles
    times_list = times.tolist()
    
    ### QUA code ###
    with program() as timeRabiProg:
    
        n = declare(int)
        t = declare(int)  # Sweeping parameter over the set of durations
        I = declare(fixed)
        Q = declare(fixed)
        I_background = declare(fixed)
        Q_background = declare(fixed)
        I_tot = declare(fixed)
        Q_tot = declare(fixed)
        
        I_stream = declare_stream()
        Q_stream = declare_stream()
        t_stream = declare_stream()
        n_stream = declare_stream()
        
        
        with for_(n, 0, n < n_avg, n + 1):
            
            save(n, n_stream)

            with for_each_(t, times_list):  # Sweep pulse duration
                
                play("gauss" * amp(amp_q_scaling), "qubit", duration=t)
                align("qubit", "rr")
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                save(I, I_stream)
                save(Q, Q_stream)
                wait(50000,"qubit")

        with stream_processing():
            I_stream.buffer(len(times)).average().save("I")
            Q_stream.buffer(len(times)).average().save("Q")
            n_stream.save('n')   
        
    I, Q = getIQ(config, timeRabiProg, showprogress=True, n_avg = n_avg)
    
    figname = f'{foldername}/figures/time_rabi_gaussamp_{round(gauss_amp * amp_q_scaling, 3)}_{len(times)}_datapoints_n_{n_avg}_' + timestring()
    datname = f'{foldername}/data/time_rabi_gaussamp_{gauss_amp * amp_q_scaling}_{len(times)}_datapoints_n_{n_avg}_' + timestring()
      
    if magphase:
        plotfunc = plotIQ_new
    
    if makeplots:
        # title = f'time rabi, gauss amp = {round(gauss_amp * amp_q_scaling, 3)}, {len(times)} datapoints, N = {n_avg} averages'
        plot_rabi(times * 4, I, xlabel = "t (ns)")
        plt.savefig(figname + "_I" + ".png")
        plt.close()
        
        plot_rabi(times * 4, Q, figname = figname + "_Q", ylabel = "Q (V)", xlabel = "t (ns)")
        plt.savefig(figname + "_Q" + ".png")
        plt.close()
        
        # plotIQ_new(4 * times, I, Q, title=title, xlabel = "t (ns)")
        # plt.savefig(figname + '.png')
        # plt.close()
        
        # plotIQ_quads(4 * times, I, Q, title=title, xlabel = "t (ns)")
        # plt.savefig(figname + '_quads.png')
        # plt.close()
       
       # save data
    with open(datname + '.csv', 'w', newline='') as f:
       writer = csv.writer(f)
       for r in [times, I, Q]:
           writer.writerow(r)
   
    
    return times, I, Q;



def powerrabi(a_min = 0.01,    # minimum amp_q scaling
               a_max = 0.5,     # maximum amp_q scaling
               da = 0.005,       # step of amp_q
               n_avg = 2000,    # number of averages
               pulse_len = 1500,  # pulse length in nanoseconds
               makeplots = True,
               foldername = foldername + "\\power_rabi"):
    
    if not os.path.isdir(foldername):
       os.mkdir(foldername);
       os.mkdir(foldername + "\\figures")
       os.mkdir(foldername + "\\data")
          
    amps = np.arange(a_min, a_max + da/2, da)
    amps_list = amps.tolist()
    pulse_len_clk = int(pulse_len/4)
    
    ### QUA code ###
    with program() as power_rabi:
    
        a = declare(fixed)
        n = declare(int)
        I = declare(fixed)
        Q = declare(fixed)
                
        I_stream = declare_stream()
        Q_stream = declare_stream()
        n_stream = declare_stream()
    
        
        with for_(n, 0, n < n_avg, n + 1):

            save(n, n_stream) 
            
            with for_each_(a, amps_list):  # Sweep pulse duration
                
                play("gauss" * amp(a), "qubit", duration = pulse_len_clk) # time in clockcycles as parameter
                align("qubit", "rr")
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                save(I, I_stream)
                save(Q, Q_stream)
                wait(50000,"qubit")

        with stream_processing():
            I_stream.buffer(len(amps)).average().save("I")
            Q_stream.buffer(len(amps)).average().save("Q")
            n_stream.save('n')
        
    I, Q = getIQ(config, power_rabi, showprogress = True, n_avg = n_avg)
    
    name = make_name(precision = 3,
                      scientific = True,
                      a_min = a_min,    # minimum amp_q scaling
                        a_max = a_max,     # maximum amp_q scaling
                        da = da,       # step of amp_q
                        n_avg = n_avg,    # number of averages
                        pulse_len = pulse_len)
    name += '_' + timestring()
    
    figname = f'{foldername}/figures/' + name
    datname = f'{foldername}/data/' + name 
       
    if makeplots:
        
        plot_rabi(amps, I, title = f"pulse length {pulse_len} ns, " )
        plt.savefig(figname + "_I.png")
        plt.close()
        
        plot_rabi(amps, Q, title = f"pulse length {pulse_len} ns, ",  ylabel = "Q (V)")
        plt.savefig(figname + "_Q.png")
        plt.close()
        # title = f'power rabi, pulse = {pulse_len} ns, N = {n_avg} averages'
        # plotIQ_new(amp_q * amps, I, Q, title=title, xlabel = "amplitude (qua units)")
        # plt.savefig(figname + ".png")
        # plt.close()
        
        # plotIQ_quads(amp_q * amps, I, Q, title=title, xlabel = "amplitude (qua units)")
        # plt.savefig(figname + '_quads.png')
        # plt.close()
        
        # save data
    with open(datname + '.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for r in [amps, I, Q]:
            writer.writerow(r)
    
        
    return amps, I, Q;


def plot_rabi(xdata, I, figname = "rabi", 
                          xlabel = "amplitude scaling", 
                          ylabel = "I (V)",
                          title = ""):
    
    p = plt.plot(xdata, I)
    
    try:
        resI = fit_sin(xdata, I);
    except RuntimeError:
        pass
    else:
        plt.plot(xdata, resI["fitfunc"](xdata))
        title += f'period = {resI["period"]}'
    finally:
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        plt.title(title)
        
        plt.tight_layout()
        
        return p



"""
"----------- time rabi ----------"

# Importing the necessary from qm
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm.qua import math
from qm import LoopbackInterface
from qm import SimulationConfig
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit
from config_MUN11 import *
from analysis_functions import convert_V_to_dBm
# import plot_functions as pf



t_min = 1
t_max = 400  # Maximum pulse duration (in clock cycles, 1 clock cycle =4 ns)
dt = 1  # timestep in clock cycles

N_max = 1000


qmm = QuantumMachinesManager()  # Reach OPX's IP address

with program() as timeRabiProg:  # Time Rabi QUA program
    I = declare(fixed)  # QUA variables declaration
    Q = declare(fixed)
    t = declare(int)  # Sweeping parameter over the set of durations
    Nrep = declare(int)  # Number of repetitions of the experiment
    I_stream = declare_stream()  # Declare streams to store I and Q components
    Q_stream = declare_stream()
    t_stream = declare_stream()
    
    with for_(Nrep, 0, Nrep < N_max, Nrep + 1):  # Do a 100 times the experiment to obtain statistics
        with for_(t, t_min, t <= t_max, t + dt):  # Sweep from 0 to 100 *4 ns the pulse duration
            
            play("gauss", "qubit", duration=t)
            align("qubit", "rr")
            measure("readout", "rr", None,
                    dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                    dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
            save(I,I_stream)
            save(Q,Q_stream)
            wait(50000,"qubit")
            # wait(50000,"rr")

    with stream_processing():
        I_stream.buffer(len(t_arr)).average().save("I")
        Q_stream.buffer(len(t_arr)).average().save("Q")
        

qm = qmm.open_qm(config)
job = qm.execute(timeRabiProg)
res_handles = job.result_handles
# res_handles.wait_for_all_values()
# a = plt.figure()
I_handle = res_handles.get("I")
Q_handle = res_handles.get("Q")
I_handle.wait_for_values(1)
Q_handle.wait_for_values(1)
while(I_handle.is_processing()):
    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    # mag = np.sqrt(I**2 + Q**2)
    # plt.figure(a)
    # plt.clf()
    # plt.plot(t_arr*4,mag)
    # plt.title('qubit spectroscopy analysis')
    # plt.xlabel("freq (GHz)")
    # plt.ylabel("Amplitude")
    # plt.pause(0.1)
    
    # plt.plot(freqs+qLO, I)
    # plt.plot(freqs+qLO, Q)
    
mag = np.sqrt(I**2 + Q**2)  
plt.figure() 
plt.plot(t_arr*4,mag)
plt.title('time-rabi')
plt.xlabel("pulse duration (ns)")
plt.ylabel("Magnitude")    

angle = np.unwrap(np.angle(I+1j*Q))
plt.figure()
plt.plot(t_arr*4,angle)
plt.title('time-rabi')
plt.xlabel("pulse duration (ns)")
plt.ylabel("Phase (rad)")




"----------- power rabi ------------"


# where figures are saved
datapath = f'D:\\Users\\lfl\\data\\MUNIN11\\res_{resFreq/1e9}GHz\\'
foldername = datapath + 'rabi_data\\'

from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from config_MUN11 import *

import matplotlib.pyplot as plt
import numpy as np

from file_utility import *
###################
# The QUA program #
###################

a_min = 0.01
a_max = 1.0
da = 0.01
amps = np.arange(a_min, a_max + da/2, da)
n_avg = 8000
detun = 0e6

with program() as power_rabi:

    n = declare(int)
    I = declare(fixed)
    I_st = declare_stream()
    Q = declare(fixed)
    Q_st = declare_stream()
    a = declare(fixed)
    update_frequency('qubit',(qbFreq-qb_LO)+detun)
    with for_(n, 0, n < n_avg, n + 1):
        with for_(a, a_min, a < a_max + da/2, a + da):
            play("gauss"*amp(a), "qubit",duration= pi_len)  # in cycles
            align("qubit", "rr")
            measure("readout", "rr", None,
                    dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                    dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
            save(I, I_st)
            save(Q, Q_st)
            wait(50000, "qubit")

    with stream_processing():
        I_st.buffer(len(amps)).average().save('I')
        Q_st.buffer(len(amps)).average().save('Q')

######################################
# Open Communication with the Server #
######################################
qmm = QuantumMachinesManager()  # Reach OPX's IP address

####################
# Simulate Program #
####################
qm = qmm.open_qm(config)

job = qm.execute(power_rabi)
res_handles = job.result_handles
# res_handles.wait_for_all_values()
# a = plt.figure()
I_handle = res_handles.get("I")
Q_handle = res_handles.get("Q")
I_handle.wait_for_values(1)
Q_handle.wait_for_values(1)

while(I_handle.is_processing()):
    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    # mag = np.sqrt(I**2 + Q**2)
    # # plt.figure()
    # plt.plot(amps, 1e3*mag,'-o')
    # # plt.plot(amps, I)
    # # plt.plot(amps, Q)
    # # plt.title('qubit spectroscopy analysis')
    # plt.xlabel("Drive Amplitude a")
    # plt.ylabel("Amplitude (mV)")
    # plt.pause(0.1)
    # plt.clf()

I = I_handle.fetch_all()
Q = Q_handle.fetch_all()    
# plt.figure()
# plt.plot(amps,I*1e3,'-o')
# plt.plot(amps,Q,'-o')
# plt.figure()
# plt.plot(amps,Q)
# plt.figure()
# plt.plot(amps,mag)

mag = np.sqrt(I**2 + Q**2)  
plt.figure() 
plt.plot(amps,mag)
plt.title('power-rabi')
plt.ylabel("Magnitude")    

angle = np.unwrap(np.angle(I+1j*Q))
plt.figure()
plt.plot(amps,angle)
plt.title('power-rabi')
plt.ylabel("Phase (rad)")


dataDict = {'amps': amps,
             'I': I,
             'Q': Q}

# datapath = filepath + datestr
# filename = 'powerRabi_detun=%sMHz_%s.csv'%(str(detun/1e6).replace('.','d'),datestr)
# filename = make_new_filename(datapath, filename, 'csv')
# dataDf = pd.DataFrame(dataDict) 

# dataDf.to_csv(datapath +'\\'+ filename)
"""
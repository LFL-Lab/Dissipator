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
    
def res_demod(I, Q):
    
    return (dual_demod.full("integW_cos", "out1", "integW_minus_sin", "out2", I),
            dual_demod.full("integW_sin", "out1", "integW_cos", "out2", Q))



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
                measure("readout", "rr", None, *res_demod(I, Q))
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
        plt.show()
        plt.close()
        
        plot_rabi(times * 4, Q, figname = figname + "_Q", ylabel = "Q (V)", xlabel = "t (ns)")
        plt.savefig(figname + "_Q" + ".png")
        plt.show()
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
                measure("readout", "rr", None, *res_demod(I, Q))
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
        plt.show()
        plt.close()
        
        plot_rabi(amps, Q, title = f"pulse length {pulse_len} ns, ",  ylabel = "Q (V)")
        plt.savefig(figname + "_Q.png")
        plt.show()
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


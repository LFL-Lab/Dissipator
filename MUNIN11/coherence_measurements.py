import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager

from config_MUN11 import *
from meas_utilities import *
from plotting import *

import sys
import os
sys.path.append("D:\Program Files\Keysight\Labber\Script")
import Labber

#########
# Connect to instrument through Labber
#########

client = Labber.connectToServer()



# SC = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='24'))
# SC.startInstrument()

###################
# The QUA program #
###################

#####################
# This program does resonator spectroscopy around a given center frequency (the estimate of the resonator
# frequency) and with a given span. This is assuming we already know roughly where we expect the resonator
# to be (e.g. using VNA), and we are looking to do a finer examination using the OPX.
#####################


res = resFreq # 7.026640e9 # estimated resonator frequency

# where to save figures
datapath = f'D:\\Users\\lfl\\data\\MUNIN11\\res_{resFreq/1e9}GHz\\'
if not os.path.isdir(datapath):
    os.mkdir(datapath);

###################
# The QUA program #
###################

def T1(n_avg = 2000, 
       t0 = 20,
       tf = 2e3, # max time in nanoseconds
       dt = 20, # time step in nanoseconds
       resettime = 200e3): # qubit reset time in nanoseconds

    args = locals() 
    foldername = makefolder(datapath, "T1")

    name = foldername + make_name(paramdict = args) + "_" + timestring()
    title = make_title(paramdict = args)

    t_min = round(t0 / 4)
    t_max = round(tf / 4)
    step = round(dt / 4)
    resettime_clk = round(resettime / 4)
    times = np.arange(t_min, t_max + step/2, step, dtype = int)
    times_list = times.tolist()
    
    with program() as T1Prog:
        
        n = declare(int)
        t = declare(int)  # Sweeping parameter over the set of durations
        I = declare(fixed)
        Q = declare(fixed)
        I_tot = declare(fixed)
        Q_tot = declare(fixed)
        
        I_stream = declare_stream()
        Q_stream = declare_stream()
        t_stream = declare_stream()
        n_stream = declare_stream()
        
        with for_(n, 0, n < n_avg, n + 1):
            
            save(n, n_stream)

            with for_each_(t, times_list):  # Sweep pulse duration
                
                with if_(t == 0):
                    
                    wait(resettime_clk, "qubit")
                    play("pi", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
                    
                with else_():
                    
                    wait(resettime_clk, "qubit")
                    play("pi", "qubit")
                    align("qubit", "rr")
                    wait(t, 'rr')
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
        

        with stream_processing():
            I_stream.buffer(len(times)).average().save("I")
            Q_stream.buffer(len(times)).average().save("Q")
            n_stream.save('n')   
        
    I, Q = getIQ(config, T1Prog, showprogress=True, n_avg = n_avg)
    
    plotIQ_quads(times * 4 / 1e3, I, Q, xlabel = "t (us)", title=title)
    plt.savefig(name + ".png")
    plt.close()
    
    datadict = {"times" : times, "I" : I, "Q" : Q}
    writedata(name + ".csv", datadict = datadict)
    
    return times * 4, I, Q
       

def ramsey(n_avg = 5000,
           t0 = 0,
           tf = 2000, # in nanoseconds
           dt = 128, # in nanoseconds
           resettime = 20e3,
           detuning = 0.0):
    
    args = locals()
    
    subfoldername = makefolder(datapath, "ramsey")
    name = subfoldername + make_name(paramdict = args) + "_" + timestring()
    title = make_title(paramdict = args)
    
    t_min = round(t0 / 4)
    t_max = round(tf / 4)
    step = round(dt / 4)
    resettime_clk = round(resettime / 4)
    times = np.arange(t_min, t_max + step/2, step, dtype = int)
    times_list = times.tolist()
    
    with program() as ramseyProg:
        
        update_frequency('qubit', (qbFreq-qb_LO) + detuning)

        n = declare(int)
        t = declare(int) 
        I = declare(fixed)
        Q = declare(fixed)
        I_tot = declare(fixed)
        Q_tot = declare(fixed)
        
        I_stream = declare_stream()
        Q_stream = declare_stream()
        n_stream = declare_stream()
        
        with for_(n, 0, n < n_avg, n + 1):
            
            save(n, n_stream)

            with for_each_(t, times_list):  # Sweep pulse duration
                
                with if_(t == 0):
                    
                    play("pi_half", "qubit")
                    play("pi_half", "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
                    wait(resettime_clk,"qubit")
                
                    
                with else_():
                    
                    play("pi_half", "qubit")
                    wait(t, "qubit")
                    # frame_rotation_2pi(phi, 'qubit')  # this was in Haimeng's code and was commented out by her,
                                                        # not sure what it's for.
                    play("pi_half", "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
                    wait(resettime_clk, "qubit")
        

        with stream_processing():
            I_stream.buffer(len(times)).average().save("I")
            Q_stream.buffer(len(times)).average().save("Q")
            n_stream.save('n')   
            
        
    I, Q = getIQ(config, ramseyProg, showprogress=True, n_avg = n_avg)
    
    plotIQ_quads(times * 4, I, Q, xlabel = "t (ns)", title=title)
    plt.savefig(name + ".png")
    plt.close()
    
    datadict = {"times" : times, "I" : I, "Q" : Q}
    writedata(name + ".csv", datadict = datadict)
    
    return datadict


def echo(n_avg = 2000,
         t0 = 0,
         tf = 6000, # in nanoseconds
         dt = 128, # in nanoseconds
         resettime = 20e3):
    
    args = locals()
    
    subfoldername = makefolder(datapath, "echo")

    name = subfoldername + make_name(paramdict = args) + "_" + timestring()
    title = make_title(paramdict = args)
    
    t_min = round(t0 / 4)
    t_max = round(tf / 4)
    step = round(dt / 4)
    resettime_clk = round(resettime / 4)
    times = np.arange(t_min, t_max + step/2, step, dtype = int)
    times_list = times.tolist()
    
    with program() as ramseyProg:
        
        n = declare(int)
        t = declare(int) 
        I = declare(fixed)
        Q = declare(fixed)
        I_tot = declare(fixed)
        Q_tot = declare(fixed)
        
        I_stream = declare_stream()
        Q_stream = declare_stream()
        n_stream = declare_stream()
        
        with for_(n, 0, n < n_avg, n + 1):
            
            save(n, n_stream)

            with for_each_(t, times_list):  # Sweep pulse duration
                
                with if_(t == 0):
                    
                    play("pi_half", "qubit")
                    play("pi", "qubit")
                    play("pi_half", "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
                    wait(resettime_clk,"qubit")
                
                    
                with else_():
                    
                    play("pi_half", "qubit")
                    wait(t/2, "qubit")
                    play("pi", "qubit")
                    wait(t/2, "qubit")
                    play("pi_half", "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
                    wait(resettime_clk, "qubit")
        

        with stream_processing():
            I_stream.buffer(len(times)).average().save("I")
            Q_stream.buffer(len(times)).average().save("Q")
            n_stream.save('n')   
            
        
    I, Q = getIQ(config, ramseyProg, showprogress=True, n_avg = n_avg)
    
    plotIQ_quads(times * 4, I, Q, xlabel = "t (ns)", title=title)
    plt.savefig(name + ".png")
    plt.close()
    
    datadict = {"times" : times, "I" : I, "Q" : Q}
    writedata(name + ".csv", datadict = datadict)
    
    return datadict
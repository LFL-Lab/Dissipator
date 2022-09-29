import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager

from config import *
from Utilities.measurement import *
from fitting import *
from plotting import *

import sys
import os
sys.path.append("D:\Program Files\Keysight\Labber\Script")
import Labber

#########
# Connect to instrument through Labber
#########

client = Labber.connectToServer()

qLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Qubit', startup = 'Get config'))
qLO.startInstrument()
qLO.setValue('Frequency', qb_LO)
qLO.setValue('Output',True)
config['elements']['qubit']['mixInputs']['lo_frequency'] = qb_LO
config['mixers']['mixer_q1'][0]['lo_frequency'] = qb_LO

rrLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Readout', startup = 'Get config'))
rrLO.startInstrument()
rrLO.setValue('Frequency', rr_LO)
rrLO.setValue('Output',True)
config['elements']['rr']['mixInputs']['lo_frequency'] = rr_LO
config['mixers']['mixer_rl1'][0]['lo_frequency'] = rr_LO



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
datapath = f'D:\\Users\\lfl\\data\\candle\\res_{resFreq/1e9}GHz\\'
if not os.path.isdir(datapath):
    os.mkdir(datapath);
    
    
def res_demod(I, Q):
    
    return (dual_demod.full("integW_cos", "out1", "integW_minus_sin", "out2", I),
            dual_demod.full("integW_sin", "out1", "integW_cos", "out2", Q))


###################
# The QUA program #
###################

def T1(n_avg = 2000, 
       t0 = 20,
       tf = 300e3, # max time in nanoseconds
       dt = 1e3, # time step in nanoseconds
       resettime = 400e3): # qubit reset time in nanoseconds

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
    
    # plotIQ_quads(times * 4 / 1e3, I, Q, xlabel = "t (us)", title=title)
    # plt.savefig(name + ".png")
    # plt.close()
    
    # datadict = {"times" : times, "I" : I, "Q" : Q}
    # writedata(name + ".csv", datadict = datadict)
    
    # make fit plots
    for (quad, data) in zip(["I", "Q"], [I, Q]):
        
        times = np.array(times) * 4 / 1e3; # times in microseconds
        data = np.array(data); 
        
        fit, fig = perform_fit(Exp, times, data, plot=True, 
                          maketitle=True,
                          title = "T1", precision = 4, 
                          xlabel = "t (us)", ylabel = quad  + " (V)")
        # fit = perform_fit(times, I, fitfunc = expCosine, plot=True, freq_0 = 0.25);
        
        fig.show()
        fig.savefig(name + "_" + quad + "fit.png")
        fig.close()
        
        print(quad + " fit:\n")
        print(fit)
        print("\n\n")
      

    
    return times, I, Q
       

def ramsey(n_avg = 2000,
           t0 = 0,
           tf = 150e3, # in nanoseconds
           dt = 250, # in nanoseconds
           resettime = 400e3,
           detuning = 20e3):
    
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
    
    for (quad, data) in zip(["I", "Q"], [I, Q]):
        
        times = np.array(times) * 4 / 1e3; # times in microseconds
        data = np.array(data); 
        
        fit, fig = perform_fit(expCosine, times, data, plot=True, 
                          maketitle=True,
                          title = "ramsey", precision = 4, 
                          xlabel = "t (us)", ylabel = quad  + " (V)")
        
        fig.show()
        fig.savefig(name + "_" + quad + "fit.png")
        fig.close()
        
        print(quad + " fit:\n")
        print(fit)
        print("\n\n")
      
    
    fit, fig = perform_fit(expCosine, times, data, plot=True, 
                      maketitle=True, freq = 0.15,
                      title = "T1", precision = 4, 
                      xlabel = "t (us)", ylabel = quad  + " (V)")
    
    
    datadict = {"times" : times, "I" : I, "Q" : Q}
    writedata(name + ".csv", datadict = datadict)
    
    return datadict


def echo(n_avg = 2000,
         t0 = 0,
         tf = 200e3, # in nanoseconds
         dt = 250, # in nanoseconds
         resettime = 400e3):
    
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
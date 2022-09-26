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

LO = client.connectToInstrument('SignalCore SC5511A Signal Generator', dict(name='10002F1D', startup = 'Get config'))
LO.startInstrument()

DA = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(interface='USB',address='26920'))
DA.startInstrument()

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
foldername = datapath + f'resonator_spec_scan_figures_' + datestring()
if not os.path.isdir(foldername):
    os.mkdir(foldername);


def run_resonator_spec(IF_min = 0.1e6, 
                       IF_max = 400e6, 
                       df = 0.1e6,
                       n_avg = 500):
    
    freqs = np.arange(IF_min, IF_max + df/2, df)
    
    ### QUA code ###
    with program() as rr_spec:
    
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)
        
        with for_(n, 0, n < n_avg, n + 1):
            
            with for_(f, IF_min, f < IF_max + df/2, f + df):
                update_frequency("rr", f)
                wait(1000, "rr")
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                save(I, I_st)
                save(Q, Q_st)

        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')
        
    I, Q = getIQ(config, rr_spec)
    
    return I, Q;

def run_resonator_spec_qubitpower(IF_min = 0.1e6, 
                                   IF_max = 400e6, 
                                   df = 0.1e6,
                                   n_avg = 500,
                                   amp_q_scaling = 1.0):
    
    freqs = np.arange(IF_min, IF_max + df/2, df)
    
    ### QUA code ###
    with program() as rr_spec:
    
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)
        
        with for_(n, 0, n < n_avg, n + 1):
            
            with for_(f, IF_min, f < IF_max + df/2, f + df):
                
                update_frequency('rr', f)
                play("const" * amp(amp_q_scaling), "qubit", duration = 5000)
                align('qubit', 'rr')
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                wait(1000, 'rr')
                save(I, I_st)
                save(Q, Q_st)

        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')
        
    I, Q = getIQ(config, rr_spec, showprogress = False, n_avg = n_avg, notify = False)
    
    return I, Q;


def scan_resonator_spec(center = res,
                        attens = np.arange(50,0,-0.5), 
                        span = 5e6, 
                        df = 0.01e6, 
                        detuning = -50e6,
                        n_avg = 100000,
                        save_traces=False):
    
    # keep the resonator LO 50 MHz above target frequency
    # ("above" implemented by minus sign on detuning)
    rrLO = center - detuning
    LO.setValue('Frequency', rrLO)
    config['elements']['rr']['mixInputs']['lo_frequency'] = rrLO
    config['mixers']['mixer_rl1'][0]['lo_frequency'] = rrLO
    
    allI, allQ = [], []
    
    for atten in tqdm(attens):
        
        # update attenuation on digital attenuation to adjust readout tone power
        DA.setValue('Attenuation', atten)
        
        # IF frequencies to be swept. These will modulate the readout LO frequency 
        # to perform a sweep from center - span/2 to center + span/2.
        IF_min = detuning - span/2
        IF_max = detuning + span/2
        
        I, Q = run_resonator_spec(IF_min, IF_max, df, n_avg)
        
        allI.append(I)
        allQ.append(Q)
        
        freqs = np.arange(IF_min, IF_max + df/2, df)
        resolution = round(df / (1e6), 3)
        if save_traces:
            p = plotIQ_new( (freqs + rrLO)/1e9, I, Q, title = f'resonator spectroscopy \n DA = {atten} dB, n_avg = {n_avg}, df = {resolution} MHz');
            plt.savefig(f'{foldername}/res_spec_DA_{atten}_dB_n_{n_avg}.png')
            plt.close()
        
    plt.close('all')
    magnitudes = mags(allI, allQ)
    hp = heatplot((freqs + rrLO)/1e9, attens, np.log10(magnitudes), xlabel = 'freq (GHz)', ylabel = 'digital atteneuation (dB)', normalize=True)
    plt.tight_layout()
    plt.savefig(f'{foldername}/res_spec_DA_{attens[1]}_to_{attens[-1]}_dB_n_{n_avg}.png')
    plt.close()    
    
    print(f"Figures saved to {foldername}")
    
    return allI, allQ, freqs + rrLO, attens



def scan_resonator_spec_qubitpower(center = res,
                                    amp_q_scalings = np.arange(0.001, 1.0, 0.01), 
                                    span = 5e6, 
                                    df = 0.01e6, 
                                    detuning = -50e6,
                                    n_avg = 100000,
                                    save_traces=False):
    
    # keep the resonator LO 50 MHz above target frequency
    # ("above" implemented by minus sign on detuning)
    rrLO = center - detuning
    LO.setValue('Frequency', rrLO)
    config['elements']['rr']['mixInputs']['lo_frequency'] = rrLO
    config['mixers']['mixer_rl1'][0]['lo_frequency'] = rrLO
    
    allI, allQ = [], []
    
    for amp_q_scaling in amp_q_scalings:
        
        # IF frequencies to be swept. These will modulate the readout LO frequency 
        # to perform a sweep from center - span/2 to center + span/2.
        IF_min = detuning - span/2
        IF_max = detuning + span/2
        
        I, Q = run_resonator_spec_qubitpower(IF_min, IF_max, df, n_avg, amp_q_scaling = amp_q_scaling)
        
        allI.append(I)
        allQ.append(Q)
        
        freqs = np.arange(IF_min, IF_max + df/2, df)
        resolution = round(df / (1e6), 3)
        if save_traces:
            p = plotIQ_new( (freqs + rrLO)/1e9, I, Q, title = f'resonator spectroscopy \n scaling = {amp_q_scaling}, n_avg = {n_avg}, df = {resolution} MHz');
            plt.savefig(f'{foldername}/res_spec_scaling_{amp_q_scaling}_n_{n_avg}.' + timestring() + 'png')
            plt.close()
        
    plt.close('all')
    magnitudes = mags(allI, allQ)
    hp = heatplot((freqs + rrLO)/1e9, amp_q_scalings, np.log10(magnitudes), xlabel = 'freq (GHz)', ylabel = 'amp_q scaling', normalize=True)
    plt.tight_layout()
    plt.savefig(f'{foldername}/res_spec_scaling_{amp_q_scalings[1]}_to_{amp_q_scalings[-1]}_n_{n_avg}_' + timestring() + '.png')
    plt.close()    
    
    print(f"Figures saved to {foldername}")
    
    return allI, allQ, freqs + rrLO, amp_q_scalings


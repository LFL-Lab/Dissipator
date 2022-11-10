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
from meas_utilities import *
from plotting import *

from slacker import sendslack



import sys
sys.path.append("D:\Program Files\Keysight\Labber\Script")
sys.path.append(r"C:\Users\lfl\measurements_test")
import Labber

from Slack_wrapper.Slack import Slack_er

slack_channel = "hpc-notifications"
sacha_id = "W018J5SNCCT"

def sendslack(slack_id = f"<@{sacha_id}> ", message = "Measurement finished."):
    Slacker = Slack_er()
    Slacker.send_message(slack_channel, slack_id + message)


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
# This module includes functions to run qubit spectroscopy over a broad range of 
# frequencies, as well as running qubit spectroscopy in a narrow range of frequencies
# while sweeping DC bias.
#####################

# where to save figures
datapath = f'D:\\Users\\lfl\\data\\candle\\res_{resFreq/1e9}GHz\\'
if not os.path.isdir(datapath):
    os.mkdir(datapath);
foldername = datapath + 'qubit_spec_figures_' + datestring()
if not os.path.isdir(foldername):
    os.mkdir(foldername);

# other parameters
saturation_dur = 25000 # time qubit saturated w/ qubit tone, in clock cycles; 12500 * 4 ns = 50 us
wait_period = 5000 # wait time between experiments, in clock cycles; 50000 * 4 ns = 200 us

res_ringdown_time = int(4e3 / 4) # ringdown time in nanoseconds / (4 ns / clockcycle)
    
"""
run_qubit_spec(
                IF_min = 0.1e6,       :: minimum IF frequency on qubit mixer 
                IF_max = 400e6,       :: maximum IF frequency on qubit mixer 
                df = 0.1e6,           :: scan resolution
                amp_q_scaling = 0.1,  :: scaling of input power to mixer
                n_avg = 1000):        :: number of averages

Sweeps IF frequencies on qubit mixer between IF_min and IF_max in order to do qubit spectroscopy.
The qubit LO frequency determines the lower bound of the sweep, i.e. if the qubit LO is 6.0 GHz,
the sweep will run from (6.0e9 + IF_min) to (6.0e9 + IF_max) in steps of df.

Returns: 
    I, Q :: lists of length = len(freqs), where freqs = np.arange(IF_min, IF_max + df/2, df),
            the averaged I and Q values corresponding to each frequency.
"""
 
def run_qubit_spec(IF_min = 0.1e6, 
                   IF_max = 400e6, 
                   df = 0.1e6,
                   amp_q_scaling = 0.1,
                   n_avg = 500,
                   notify = False):
    
    freqs = np.arange(IF_min, IF_max + df/2, df)
    
    ### QUA code ###
    with program() as qubit_spec:
    
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)
        
        with for_(n, 0, n < n_avg, n + 1):

            with for_(f, IF_min, f < IF_max + df/2, f + df):
                update_frequency("qubit", f) # updates IF frequency on qubit mixer
                play("const" * amp(amp_q_scaling), "qubit", duration = saturation_dur)
                align("qubit", "rr")
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                
                save(I, I_st)
                save(Q, Q_st)
                wait(wait_period, "rr")
    
        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')
        
    I, Q = getIQ(config, qubit_spec)
    
    return I, Q, freqs;

def setLO(value, element = "qubit"):
    
    if element == "qubit":
    
        qLO.setValue('RF1 frequency', value)
        config['elements']['qubit']['mixInputs']['lo_frequency'] = value
        config['mixers']['mixer_q1'][0]['lo_frequency'] = value
     
    elif element == "rr":
        
        rrLO.setValue('Frequency', value)
        config['elements']['rr']['mixInputs']['lo_frequency'] = value
        config['mixers']['mixer_rl1'][0]['lo_frequency'] = value
        
    else:
        raise Exception("Please choose 'qubit' or 'rr' as the element.")


"modify this for measure, pulse, measure -- 'on-off measurement / background subtraction / lockin measurement' "
def run_qubit_spec_onoff(IF_min = 0.1e6, 
                   IF_max = 200e6, 
                   df = 1.0e6,
                   amp_q_scaling = 1.0,
                   n_avg = 500,
                   notify = False,
                   saveplot = False):
    
    freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
    freqs_list = freqs.tolist()
    
    ### QUA code ###
    with program() as qubit_spec:
    
        n = declare(int)
        f = declare(int)
        I = declare(fixed)
        Q = declare(fixed)
        I_background = declare(fixed)
        Q_background = declare(fixed)
        I_tot = declare(fixed)
        Q_tot = declare(fixed)
        
        I_stream = declare_stream()
        Q_stream = declare_stream()
        n_stream = declare_stream()

        
        with for_(n, 0, n < n_avg, n + 1):
            
            save(n, n_stream)

            with for_each_(f, freqs_list):
                update_frequency("qubit", f) # updates IF frequency on qubit mixer
                
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I_background),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q_background))
                wait(res_ringdown_time, "rr")
                align("rr", "qubit") # wait for operations on resonator to finish before playing qubit pulse

                
                play("const" * amp(amp_q_scaling), "qubit", duration = saturation_dur)
                wait(50, "qubit") # WHY is this here? old testing code...
                align("qubit", "rr") # wait for operations on qubit to finish before playing readout pulse
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                
                assign(I_tot, I - I_background)
                assign(Q_tot, Q - Q_background)
                save(I_tot, I_stream)
                save(Q_tot, Q_stream)
                wait(wait_period, "rr")
                
         
        with stream_processing():
            I_stream.buffer(len(freqs)).average().save('I')
            Q_stream.buffer(len(freqs)).average().save('Q')
            n_stream.save('n')    
            
    I, Q = getIQ(config, qubit_spec, showprogress=True, n_avg = n_avg, notify = notify)
    
    if saveplot:
        name = f'{foldername}/'
        LOfreq = qb_LO #qLO.getValue('RF1 frequency')
        name += make_name(precision = 3,
                          scientific = True,
                          minfreq = IF_min + LOfreq,
                          maxfreq = IF_max + LOfreq,
                          df = df, 
                          amp_q = amp_q * amp_q_scaling,
                          n = n_avg,
                          LO = LOfreq)
        
        name += "_upper_" + timestring()
        
        plt.clf()
        p = plotIQ_both(I, Q, freqs + LOfreq, title = f"qubit spectroscopy, amp_q = {amp_q * amp_q_scaling}");          
        plt.savefig(name + '.png')
        plt.close()
        
    
    return I, Q, freqs;


def plot_qubit_spec(freqs, I, Q, LOfreq, IF_min, IF_max, amp_q_scaling, n_avg, subtract_background = False):
    
    rminfreq = round((LOfreq + IF_min)/1e9, 4)
    rmaxfreq = round((LOfreq + IF_max)/1e9, 4)
    ramp_q = round(amp_q_scaling * amp_q, 5)
    name = f'{foldername}/qubit_spec_{rminfreq}_to_{rmaxfreq}_GHz_amp_q_{ramp_q}_n_{n_avg}'
    if subtract_background:
        name += "subtract_background"

    # create intermediate plots
    p = plotIQ_both(I, Q, freqs + LOfreq, title = "qubit spectroscopy");          
    plt.savefig(name + '.png')
    plt.close()
    
    # save intermediate data
    with open(name + '.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for r in [freqs + LOfreq, I, Q]:
            writer.writerow(r)

def run_qubit_spec_onoff_lower(IF_min = -200e6, 
                   IF_max = -325e6, 
                   df = -.25e6,
                   amp_q_scaling = .85,
                   n_avg = 2000,
                   notify = False,
                   saveplot=True):
    
    freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
    freqs_list = freqs.tolist()
    
    ### QUA code ###
    with program() as qubit_spec:
    
        n = declare(int)
        f = declare(int)
        I = declare(fixed)
        Q = declare(fixed)
        I_background = declare(fixed)
        Q_background = declare(fixed)
        I_tot = declare(fixed)
        Q_tot = declare(fixed)
        
        I_stream = declare_stream()
        Q_stream = declare_stream()
        n_stream = declare_stream()
        
        with for_(n, 0, n < n_avg, n + 1):
            
            save(n, n_stream)

            with for_each_(f, freqs_list):
                update_frequency("qubit", f) # updates IF frequency on qubit mixer
                
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I_background),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q_background))
                wait(res_ringdown_time, "rr")
                align("rr", "qubit") # wait for operations on resonator to finish before playing qubit pulse

                
                play("const" * amp(amp_q_scaling), "qubit", duration = saturation_dur)
                align("qubit", "rr") # wait for operations on qubit to finish before playing readout pulse
                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                
                assign(I_tot, I - I_background)
                assign(Q_tot, Q - Q_background)
                save(I_tot, I_stream)
                save(Q_tot, Q_stream)
                wait(wait_period, "rr")
    
        with stream_processing():
            I_stream.buffer(len(freqs)).average().save('I')
            Q_stream.buffer(len(freqs)).average().save('Q')
            n_stream.save('n')            

    I, Q = getIQ(config, qubit_spec, showprogress=True, n_avg = n_avg, notify = notify)
    
    if saveplot:
        name = f'{foldername}/'
        LOfreq = qb_LO #qLO.getValue('RF1 frequency')
        name += make_name(precision = 3,
                          scientific = True,
                          minfreq = IF_min + LOfreq,
                          maxfreq = IF_max + LOfreq,
                          df = df, 
                          amp_q = amp_q * amp_q_scaling,
                          n = n_avg,
                          LO = LOfreq)
        name += "_lower_" + timestring()
        
        plt.clf()
        p = plotIQ_both(I, Q, freqs + LOfreq, title = f"qubit spectroscopy, amp_q = {amp_q * amp_q_scaling}");          
        plt.savefig(name + '.png')
    
    return I, Q, freqs;

    

"""
scan_qubit_spec(
        minfreq = 3e9,          :: minimum frequency of sweep
        maxfreq = 6e9,          :: maximum frequency of sweep
        chunksize = 350e6,      :: size of IF 'chunk'; LO will be stepped by chunksize
        detuning = 50e6,        :: minimum detuning of qubit LO from target frequency
        df = 0.1e6,             :: resolution of sweep
        amp_q_scaling = 1.0,    :: power scaling on qubit tone
        n_avg = 1000):          :: number of averages

Performs qubit spectroscopy over a wide range between minfreq and maxfreq.
Sweeps qubit LO in steps of chunksize, then varies the qubit IF to do finer sweeps 
in each of these steps. 

Generates plots of magnitude and phase of the measured signal, as functions of frequency, 
and saves to foldername. Intermediate plots are generated for each qubit LO value, and a final plot
displaying all the data is output at the end.

Returns: nothing

NOTE: Needs modifications for lower sideband sweeps
"""
    
def scan_qubit_spec(minfreq = 2.8e9, 
                    maxfreq = 6.20e9,   
                    chunksize = 200e6,
                    detuning = 10e6, 
                    df = 2.5e6,
                    amp_q_scaling = 1.0,
                    n_avg = 2000,
                    foldername = foldername,
                    subtract_background = True,
                    lower_sideband = False,
                    notify = False):
    
    if subtract_background and lower_sideband:
        run = run_qubit_spec_onoff_lower # may need further modifications to step through LO properly
    elif subtract_background and not lower_sideband:
        run = run_qubit_spec_onoff
    elif not subtract_background and not lower_sideband:
        run = run_qubit_spec
    else:
        raise Exception('Lower sideband without background subtraction not yet implemented.')
    
    # make folder for saving data
    if not os.path.isdir(foldername):
        os.mkdir(foldername);
        
    if not os.path.isdir(foldername + '/intermediate_data'):
        os.mkdir(foldername + '/intermediate_data')
    
    if chunksize + detuning > 400e6:
        raise Exception('The maximum IF is 400 MHz. Please choose a chunksize and detuning such that'
                        'chunksize + detuning <= 400e6.')

    # fill these with I, Q values from each step of the sweep
    allI, allQ = [], []
    
    qLO_range = np.arange(minfreq - detuning, maxfreq - detuning, chunksize)
    
    for i, LOfreq in enumerate(qLO_range):
        print(f"The current LO value is {LOfreq/1e9} GHz. ({i} / {len(qLO_range)})\n")
        
        # update frequency output by Signal Generator, and update QM configuration
        setLO(LOfreq, element = "qubit")
        
        # IF frequencies to be swept. These modulate the LO frequency to perform a sweep from
        # qb_LO + IF_min to qb_LO + IF_max, with resolution df.
        IF_min = detuning + df
        IF_max = min(detuning + chunksize, maxfreq - LOfreq) # stay within the bounds of the intended sweep
        freqs = np.arange(IF_min, IF_max + df/2, df)
        
        
        # run qubit spectroscopy in the specified range
        I, Q, freqs = run(IF_min, IF_max, df, amp_q_scaling, n_avg, notify = notify)
        
        # concatenate new values of freqs, I, and Q to lists
        allI = np.concatenate((allI, I))
        allQ = np.concatenate((allQ, Q))
        
        # file naming conventions
        rminfreq = round((LOfreq + IF_min)/1e9, 4)
        rmaxfreq = round((LOfreq + IF_max)/1e9, 4)
        ramp_q = round(amp_q_scaling * amp_q, 5)
        name = f'{foldername}/intermediate_data/qubit_spec_{rminfreq}_to_{rmaxfreq}_GHz_amp_q_{ramp_q}_n_{n_avg}'
        if subtract_background:
            name += "subtract_background"
   
        # create intermediate plots
        p = plotIQ_both(I, Q, freqs + LOfreq, title = "qubit spectroscopy");          
        plt.savefig(name + '.png')
        plt.close()
        
        # save intermediate data
        with open(name + '.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            for r in [freqs + LOfreq, I, Q]:
                writer.writerow(r)
            
    
    # file naming conventions
    rminfreq = round(minfreq/1e9, 4)
    rmaxfreq = round(maxfreq/1e9, 4)
    ramp_q = round(amp_q_scaling * amp_q, 5)
    name = f'{foldername}/qubit_spec_{rminfreq}_to_{rmaxfreq}_GHz_amp_q_{ramp_q}_n_{n_avg}_df_{df/1e6}MHz'
    if subtract_background:
        name += "subtract_background"
    
    # generate plot of all frequencies
    allfreqs = np.arange(minfreq, maxfreq, df)
    p = plotIQ_both(allI, allQ, allfreqs, title = f'qubit spectroscopy df = {df/1e6} MHz, amp_q = {ramp_q}, n = {n_avg}')
    plt.savefig(name + '.png')
    plt.close()
    
    # save data from all frequencies
    with open(name + '.csv.', 'w', newline='') as f:
        writer = csv.writer(f)
        for r in [allfreqs, allI, allQ]:
            writer.writerow(r)
            
    if notify:
        sendslack(message = ":tada: Finished running spectroscopy experiment.")



def import_qubit_spec_scan(csvfilename):
    
    with open(csvfilename, 'r') as f:
        
        reader = csv.reader(f)
        freqs = np.array([float(i) for i in next(reader)])
        I = np.array([float(i) for i in next(reader)])
        Q = np.array([float(i) for i in next(reader)])
        
    return freqs, I, Q


"""
New version of scan_qubit_spec_FF that can scan large ranges of frequencies.
Can delete previous version once we confirm this is working.
"""
def sweep_FF(currents = np.arange(0.0, 1000e-6, 100e-6),
             minfreq = 4.01e9, 
             maxfreq = 6.20e9,
             df = 1e6,
             amp_q_scaling = 1.0,
             n_avg = 1000,
             foldername = foldername + '\\FF_sweep2',
             subtract_background = True,
             lower_sideband = False):
    
        if max(currents) > 3.01e-3:
            raise Exception("Maximum current setting is greater than 3.01 mA.")
    
    
        # make folder for saving data
        if not os.path.isdir(foldername):
            os.mkdir(foldername);
            
        # start DC source
        SC = client.connectToInstrument('Keithley 2400 SourceMeter', dict(interface='GPIB',address='24'))
        SC.startInstrument()
        SC.setValue('Output on', True)

        
        for current in currents:
            
            current_in_uA = round(current * 1e6, 5)
            print(f'\nThe intended current on the Keithley is {current_in_uA} uA.\n')
            
            # set the Keithley current value
            SC.setValue('Source current', current)
            voltagevalue = round(SC.getValue('Voltage') * 1e3, 5)
            print(f'\nThe measured voltage on the Keithley is {voltagevalue} mV\n')

            
            # run sweep at this current value; scan_qubit_spec automatically saves data to csv
            subfolder = f'DC_bias_{current_in_uA}_uA'
            scan_qubit_spec(minfreq = minfreq, maxfreq = maxfreq, df = df, amp_q_scaling = amp_q_scaling, n_avg = n_avg, foldername = foldername + '/' + subfolder, subtract_background = subtract_background, lower_sideband = lower_sideband)
            
            
"""            
Plots sweep data generated by sweep_FF. Takes as argument foldername, the same as input to sweep_FF.
"""            
def plot_sweep_data(foldername = foldername + "\\FF_sweep2", 
                    subfoldername = lambda val: f'DC_bias_{val}_uA',
                    sweepvals = np.arange(-2000.0, -990.0, 100.0),
                    **kwargs):
    
    freqs, sweepvals, maglist, phaselist = import_sweep_data(foldername, sweepvals = sweepvals);
    
    hp = heatplot(sweepvals, freqs/1e9, np.transpose(maglist), xlabel = "current (uA)", ylabel = "freq (GHz)", colorbarlabel = "magnitude", **kwargs)
    plt.savefig(f'{foldername}\\magplot.png')
    plt.close()
    
    
def import_sweep_data(foldername, subfoldername = lambda val: f'DC_bias_{val}_uA',
                           sweepvals = np.arange(0.0, 3000.0, 100.0)):
    
    maglist = []
    phaselist = []
    
    for val in sweepvals:
        
        subfolder = subfoldername(val)
        path = foldername + "\\" + subfolder
        files = os.listdir(path)
        
        for filename in files:
            
            if '.csv' in filename:
                csvfile = filename
                break
        
        freqs, I, Q = import_qubit_spec_scan(path + '/' + csvfile)
        maglist.append(mag(I, Q))
        phaselist.append(phase(I, Q))
        
        
    return freqs, sweepvals, np.array(maglist), np.array(phaselist)
        
    
    
        

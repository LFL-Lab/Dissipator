"""
Created on Mon Oct 4 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
"""

from tqdm import tqdm
from qm import generate_qua_script
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from slacker import sendslack
# from Utilities.data import *
from config import *
import plot_functions as pf
import os
from datetime import datetime
from meas_utilities import res_demod
import csv
import glob
import time
import instrument_init as inst
# from Resonator import Resonator as res
import numpy as np

electrical_delay = config['elements']['rr']['time_of_flight']
#%% play_pulses
def play_pulses():

    with program() as play_pulses:
        with infinite_loop_():
            play("const", 'qubit',duration=100)
            play("readout", "rr", duration=100)

    qmm = QuantumMachinesManager()
    qm = qmm.open_qm(config)
    job = qm.execute(play_pulses)

    return qm

#%% punchout
def punchout(df = 0.1e6,
             span = 20e6,
             n_avg = 500,
             atten_range = [10,30],
             atten_step = 0.1,
             res_freq = [6,7.2],
             res_ringdown_time = int(4e3)):
    """
    Executes punchout measurement for list of resonator frequencies

    Args:
        df (TYPE, optional): DESCRIPTION. Defaults to 0.1e6.
        n_avg (TYPE, optional): DESCRIPTION. Defaults to 500.
        atten_range (TYPE, optional): [min,max] attenuation values.
        res_freq (TYPE, optional): list with all the resonator frequencies in GHz. Defaults to [6e9,7.2e9].
        res_ringdown_time (TYPE, optional): DESCRIPTION. Defaults to int(4e3).

    Returns:
        None.

    """
    try:
        list_of_files = glob.glob(r'D:\weak_measurements\spectroscopy\\resonator_spec\punchout\*.csv')
        latest_file = max(list_of_files, key=os.path.getctime)
        iteration = int(latest_file[-7:-4].lstrip('0')) + 1
    except:
        iteration = 1

    attenuation_range = np.arange(atten_range[0],atten_range[1],step=atten_step)
    freq_arr = np.zeros((int(span/df),len(res_freq)))
    I = np.zeros((len(attenuation_range),int(span/df),len(res_freq)))
    Q = np.zeros((len(attenuation_range),int(span/df),len(res_freq)))
    mag = np.zeros((len(attenuation_range),int(span/df),len(res_freq)))
    magData = np.zeros((len(attenuation_range),int(span/df)))

    i = 0
    for fc in res_freq:
        print(f'Measuring resonator at {fc*1e-9} GHz')
        f_LO = fc - span/2
        j = 0
        for a in attenuation_range:
            print(f'Attenuation = {a} dB')
            dataI,dataQ,freqs,job = resonator_spec(f_LO=f_LO,atten=a,IF_min=df,IF_max=span,df=df,n_avg=n_avg,res_ringdown_time=res_ringdown_time,savedata=False,fit=False)
            freq_arr[:,i] = freqs
            I[j,:,i] = dataI
            Q[j,:,i] = dataQ
            mag[j,:,i] = np.abs(dataI+1j*dataQ)
            j += 1
        chi= freqs[np.argmin(mag[0,:,i])] - freqs[np.argmin(mag[-1,:,i])]
        print(f'Dispersive shift for resonator at {round(fc*1e-9,5)} GHz: {round(0.5*chi/np.pi*1e-3,1)} kHz')
        pf.heatplot(xdata=np.around(freqs*1e-9,4),ydata=attenuation_range,data=pf.Volt2dBm(mag[:,:,i]),xlabel='Frequency (GHz)',ylabel='Attenuation (dB)',cbar_label='Magnitude (dBm)')
        i += 1

    exp_dict = {'date/time':    datetime.now(),
               'nAverages': n_avg,
                     'w_LO': pars['rr_LO'],
                     'attenuation': attenuation_range,
            'wait_period':  res_ringdown_time,
            }

    # save data
    with open(f"D:\weak_measurements\spectroscopy\\resonator_spec\data_{iteration:03d}.csv","w") as datafile:
        writer = csv.writer(datafile)
        writer.writerow(exp_dict.keys())
        writer.writerow(exp_dict.values())
        writer.writerow(freq_arr)
        writer.writerow(I)
        writer.writerow(Q)

    return I, Q, freq_arr, job


#%% run_scan
def run_scan(df = 0.1e6,
             n_avg = 500,
             element='resonator',
             chunksize = 200e6,
             attenuation=20,
             lo_min = 6e9,
             lo_max = 7e9,
             amp_q_scaling = 1,
             saturation_dur = 20e3,
             showprogress=False,
             res_ringdown_time = int(4e3)):
    """
    Scans a broad range of frequencies in search for qubits/resonators

    Args:
        IF_min (TYPE, optional): DESCRIPTION. Defaults to 0.1e6.
        IF_max (TYPE, optional): DESCRIPTION. Defaults to 400e6.
        df (TYPE, optional): DESCRIPTION. Defaults to 0.1e6.
        n_avg (TYPE, optional): DESCRIPTION. Defaults to 500.
        res_ringdown_time (TYPE, optional): DESCRIPTION. Defaults to int(4e3).

    Returns:
        I (TYPE): DESCRIPTION.
        Q (TYPE): DESCRIPTION.
        freq_arr: DESCRIPTION.
        TYPE: DESCRIPTION.

    """

    try:
        list_of_files = glob.glob(r'D:\weak_measurements\spectroscopy\{element}_spec\*.csv')
        latest_file = max(list_of_files, key=os.path.getctime)
        iteration = int(latest_file[-7:-4].lstrip('0')) + 1
    except:
        iteration = 1

    freq_arr = []
    I = []
    Q = []

    if lo_min != lo_max:
        numchunks = int((lo_max-lo_min)/chunksize)
        lo_list = [i*chunksize+lo_min for i in range(numchunks)]
    else:
        numchunks = 1
        lo_list = [lo_min]

    for f in lo_list:
        if element == 'resonator':
            dataI,dataQ,freqs,job = resonator_spec(f_LO=f,atten=attenuation,IF_min=df,IF_max=chunksize,df=df,n_avg=n_avg,res_ringdown_time=res_ringdown_time,savedata=False)
        elif element == 'qubit':
            dataI,dataQ,freqs,job = qubit_spec(f_LO=f,amp_q_scaling=amp_q_scaling,saturation_dur=saturation_dur,atten=attenuation,IF_min=df,IF_max=chunksize,df=df,n_avg=n_avg,res_ringdown_time=res_ringdown_time,showprogress=showprogress,savedata=False)
        freq_arr.extend(freqs)
        I.extend(dataI)
        Q.extend(dataQ)

    exp_dict = {'date/time':    datetime.now(),
               'nAverages': n_avg,
                     'w_LO': pars['rr_LO'],
                     'attenuation': attenuation,
            'wait_period':  res_ringdown_time,
            }

    # save data
    with open(f"D:\weak_measurements\spectroscopy\\resonator_spec\data_{iteration:03d}.csv","w") as datafile:
        writer = csv.writer(datafile)
        writer.writerow(exp_dict.keys())
        writer.writerow(exp_dict.values())
        writer.writerow(freq_arr)
        writer.writerow(I)
        writer.writerow(Q)

    return I, Q, freq_arr, job

#%% resonator_spec
def resonator_spec(IF_min = 0.1e6,
                   f_LO = 7e9,
                   IF_max = 400e6,
                   df = 0.1e6,
                   atten = 10,
                   n_avg = 500,
                   res_ringdown_time = int(4e3),
                   port_type = 'notch',
                   fit=True,
                   savedata=True):
    """


    Args:
        IF_min (TYPE, optional): DESCRIPTION. Defaults to 0.1e6.
        f_LO (TYPE, optional): DESCRIPTION. Defaults to 7e9.
        IF_max (TYPE, optional): DESCRIPTION. Defaults to 400e6.
        df (TYPE, optional): DESCRIPTION. Defaults to 0.1e6.
        atten (TYPE, optional): DESCRIPTION. Defaults to 10.
        n_avg (TYPE, optional): DESCRIPTION. Defaults to 500.
        res_ringdown_time (TYPE, optional): DESCRIPTION. Defaults to int(4e3).
        port_type (TYPE, optional): DESCRIPTION. Defaults to 'notch'.
        fit (TYPE, optional): DESCRIPTION. Defaults to True.
        plot (TYPE, optional): DESCRIPTION. Defaults to True.
        savedata (TYPE, optional): DESCRIPTION. Defaults to True.

    Returns:
        I (TYPE): DESCRIPTION.
        Q (TYPE): DESCRIPTION.
        TYPE: DESCRIPTION.
        TYPE: DESCRIPTION.

    """

    try:
        list_of_files = glob.glob(r'D:\weak_measurements\spectroscopy\resonator_spec\*.csv')
        latest_file = max(list_of_files, key=os.path.getctime)
        iteration = int(latest_file[-7:-4].lstrip('0')) + 1
    except:
        iteration = 1

    freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
    # freqs_list = freqs.tolist()
    # set attenuation and change rr_LO freq
    inst.set_attenuator(attenuation=atten)
    pars['rr_LO'] = f_LO
    with open("pars.json", "w") as outfile:
        json.dump(pars, outfile)
    inst.set_rr_LO(f_LO)
    ### QUA code ###
    with program() as rr_spec:

        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)
        n_stream = declare_stream()

        with for_(n, 0, n < n_avg, n + 1):

            # with for_each_(f, freqs_list):
            with for_(f, IF_min, f < IF_max + df/2, f + df):

                update_frequency("rr", f)
                wait(res_ringdown_time, "rr")
                measure("readout", "rr", None,*res_demod(I, Q))
                save(I, I_st)
                save(Q, Q_st)

            save(n,n_stream)

        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')
            n_stream.save('n')

    datadict,job = get_results(config, rr_spec,result_names=["I","Q","n"],showprogress=False)

    I = datadict["I"]
    Q = datadict["Q"]
    freq_arr = freqs+pars['rr_LO']

    if fit:
        fc,fwhm = pf.fit_res(freq_arr,np.abs(I+1j*Q))
        pf.spec_plot(freq_arr,I,Q,attenuation=atten,df=df,iteration=iteration,element='resonator',fwhm=fwhm,fc=fc)
        print(f'Resonant Frequency: {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-6} MHz\nkappa = {2*np.pi*fwhm*1e-6:.3f} MHz')

    exp_dict = {'date/time':    datetime.now(),
               'nAverages': n_avg,
                     'w_LO': pars['rr_LO'],
            'wait_period':  res_ringdown_time,
            }
    if savedata:
        # save data
        with open(f"D:\weak_measurements\spectroscopy\\resonator_spec\data_{iteration:03d}.csv","w") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(exp_dict.keys())
            writer.writerow(exp_dict.values())
            writer.writerow(freqs)
            writer.writerow(I)
            writer.writerow(Q)

    return I, Q, freqs+pars['rr_LO'], job;

#%% qubit_spec
def qubit_spec(f_LO = 5e9,
                   IF_min = 0.1e6,          # min IF frequency
                   IF_max = 400e6,          # max IF frequency
                   df = 0.1e6,              # IF frequency step
                   rr_freq = 6e9,       #resonator frequency
                   amp_q_scaling = 0.1,     # prefactor to scale default "const" qubit tone, amp_q
                   n_avg = 500, # number of averages
                   atten = 10, # readout attenuation
                   saturation_dur = int(20e3),   # time qubit saturated w/ qubit tone, in ns
                   wait_period = int(40e3),      # wait time between experiments, in ns
                   res_ringdown_time = int(4e3), # resonator ringdown time for onoff measurement
                   on_off =  True,          # background subtraction
                   notify = False,
                   showprogress=False,
                   savedata=True):         # create notification on Slack when measurement finishes

    # create list of frequencies for QUA to pull from
    # NOTE: the for_each_ loop can create latency, according to QM, so we may
    # want to replace it eventually with a qualang_tools/loop functionality
    # which pulls from a python array
    #
    # SEE: https://github.com/qua-platform/py-qua-tools/tree/main/qualang_tools/loops
    try:
        list_of_files = glob.glob(r'D:\weak_measurements\spectroscopy\qubit_spec\*.csv')
        latest_file = max(list_of_files, key=os.path.getctime)
        iteration = int(latest_file[-7:-4].lstrip('0')) + 1
    except:
        iteration = 1

    freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
    saturation_dur = int(saturation_dur)
    inst.set_attenuator(attenuation=atten)

    with open('pars.json', 'r') as openfile:
        pars_dict = json.load(openfile)

    pars['qb_LO'] = f_LO
    pars['rr_LO'] = pars['rr_freq'] - 50e6

    with open("pars.json", "w") as outfile:
        json.dump(pars, outfile)
    inst.set_qb_LO(f_LO)
    inst.set_rr_LO(pars['rr_LO'])
    ### QUA code ###
    # program() delimits the QUA code, which has its own syntax separate from python
    # This code defines a QM job 'QubitSpecProg', which is later run in 'getIQ'
    with program() as QubitSpecProg:

        # first, we declare special QUA types
        n = declare(int) # averaging iterable
        f = declare(int) # frequency iterable
        I = declare(fixed)
        Q = declare(fixed)
        I_stream = declare_stream()
        Q_stream = declare_stream()
        n_stream = declare_stream()
        if on_off:
            I_background = declare(fixed)
            Q_background = declare(fixed)
            I_tot = declare(fixed)
            Q_tot = declare(fixed)
        # loop over n_avg iterations
        with for_(n, 0, n < n_avg, n + 1):

            # loop over list of IF frequencies
            with for_(f, IF_min, f < IF_max + df/2, f + df):
                # update IF frequency going into qubit mixer
                update_frequency("qubit", f)
                # measure background
                if on_off:
                    measure("readout", "rr", None, *res_demod(I_background, Q_background))
                    wait(res_ringdown_time, "rr")
                    align("rr", "qubit") # wait for operations on resonator to finish before playing qubit pulse
                # play qubit pulse and measure
                play("const" * amp(amp_q_scaling), "qubit", duration = saturation_dur)
                align("qubit", "rr") # wait for operations on resonator to finish before playing qubit pulse
                measure("readout", "rr", None, *res_demod(I, Q))
                # subtract background and save to stream
                if on_off:
                    assign(I_tot, I - I_background)
                    assign(Q_tot, Q - Q_background)
                    save(I_tot, I_stream)
                    save(Q_tot, Q_stream)
                else:
                    save(I, I_stream)
                    save(Q, Q_stream)
                # wait some time before continuing to next IF frequency
                wait(wait_period, "rr")
            save(n, n_stream)

        # average data over iterations and save to stream
        with stream_processing():
            I_stream.buffer(len(freqs)).average().save('I')
            Q_stream.buffer(len(freqs)).average().save('Q')
            n_stream.save('N')

    # execute 'QubitSpecProg' using configuration settings in 'config'
    # fetch averaged I and Q values that were saved
    datadict, job = get_results(config, QubitSpecProg, result_names = ["I", "Q", "N"], nPoints=n_avg,showprogress=showprogress, notify = notify)
    I = datadict["I"]
    Q = datadict["Q"]
    freq_arr = freqs+pars['qb_LO']

    pf.spec_plot(freq_arr,I,Q,iteration=iteration,element='qubit',find_peaks=True)
    # print(f'Qubit Frequency: {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-6} MHz\nkappa = {2*np.pi*fwhm*1e-6:.3f} MHz')

    if savedata:
        exp_dict = {'date/time':     datetime.now(),
                   'nAverages': n_avg,
                         'A_d':     amp_q_scaling,
                         'w_LO': pars['qb_LO'],
                'wait_period':  wait_period,
                }
        # save data
        with open(f"D:\weak_measurements\spectroscopy\qubit_spec\data_{iteration:03d}.csv","w") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(exp_dict.keys())
            writer.writerow(exp_dict.values())
            writer.writerow(freqs)
            writer.writerow(I)
            writer.writerow(Q)

    return I, Q, freqs, job;

#%% pulse_exp
def pulse_exp(exp,
       n_avg = 2000,
       t0 = 0,         # minimum pulse duration in nanoseconds
       tf = 10e3,    # maximum pulse duration in nanoseconds
       dt = 500,        # step of sweep in nanoseconds
       amp_q_scaling = 1,
       fit = True,
       plot = True,
       detuning = 0e6,
       resettime = 400e3):
    """

    Args:
        exp (str):   What experiment to run. Options are: 'rabi','ramsey','T1','echo'.
        n_avg (TYPE, optional): DESCRIPTION. Defaults to 2000.
        t0 (TYPE, optional): minimum pulse duration in nanoseconds. Defaults to 0.
        tf (TYPE, optional): maximum pulse duration in nanoseconds. Defaults to 10e3.
        dt (TYPE, optional): sequence step size in nanoseconds. Defaults to 500.
        amp_q_scaling (TYPE, optional): DESCRIPTION. Defaults to 1.
        plot (TYPE, optional): Whether to plot and fit the data. Defaults to True.
        detuning (float): detuning from fL0-fIF in Hz.
        resettime (TYPE, optional): waiting time between experiments. Defaults to 400e3.

    Returns:
        times (TYPE): DESCRIPTION.
        I (TYPE): DESCRIPTION.
        Q (TYPE): DESCRIPTION.

    """

    try:
        list_of_files = glob.glob(f'D:\weak_measurements\{exp}\*.csv')
        latest_file = max(list_of_files, key=os.path.getctime)
        iteration = int(latest_file[-7:-4].lstrip('0')) + 1
    except:
        iteration = 1


    t_min = round(t0 / 4)
    t_max = round(tf / 4)
    step = round(dt / 4)
    resettime_clk = round(resettime / 4)
    t_arr = np.arange(t_min, t_max + step/2, step, dtype = int)
    times_list = t_arr.tolist()


    if exp == 'rabi':
        with program() as prog:

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
                    wait(resettime_clk,"qubit")

            with stream_processing():
                I_stream.buffer(len(t_arr)).average().save("I")
                Q_stream.buffer(len(t_arr)).average().save("Q")
                n_stream.save('n')

    elif exp == 'ramsey':

        with program() as prog:

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
                       measure("readout", "rr", None,*res_mod(I,Q))
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
                       measure("readout", "rr", None, *res_demod(I, Q))
                       save(I, I_stream)
                       save(Q, Q_stream)
                       wait(resettime_clk, "qubit")


           with stream_processing():
               I_stream.buffer(len(times)).average().save("I")
               Q_stream.buffer(len(times)).average().save("Q")
               n_stream.save('n')

    elif exp == 'echo':
        with program() as prog:

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
                        measure("readout", "rr", None,*res_mod(I,Q))
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
                        measure("readout", "rr", None, *res_demod(I, Q))
                        save(I, I_stream)
                        save(Q, Q_stream)
                        wait(resettime_clk, "qubit")


            with stream_processing():
                I_stream.buffer(len(times)).average().save("I")
                Q_stream.buffer(len(times)).average().save("Q")
                n_stream.save('n')

    elif exp == 'T1':
        with program() as prog:

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
                        measure("readout", "rr", None, *res_mod(I,Q))
                        save(I, I_stream)
                        save(Q, Q_stream)

                    with else_():

                        wait(resettime_clk, "qubit")
                        play("pi", "qubit")
                        align("qubit", "rr")
                        wait(t, 'rr')
                        measure("readout", "rr", None,*res_mod(I,Q))
                        save(I, I_stream)
                        save(Q, Q_stream)


            with stream_processing():
                I_stream.buffer(len(times)).average().save("I")
                Q_stream.buffer(len(times)).average().save("Q")
                n_stream.save('n')

    I, Q,job = getIQ(config, prog, showprogress=True, n_avg = n_avg)

    t_arr = np.array(t_arr) * 4 / 1e3; # times in microseconds
    ydata = abs(I+1j*Q)

    if plot:
        fitted_pars, error = pf.fit_data(t_arr,ydata,sequence=exp,dt=t_arr[-1]*1e-6/len(t_arr))
        pf.plot_data(t_arr,ydata,sequence=exp,fitted_pars=fitted_pars,nAverages=n_avg,
                     qubitDriveFreq=qb_LO-qb_IF,amplitude_hd=amp_q_scaling*0.45,iteration=iteration)

    exp_dict = {'date/time':     datetime.now(),
               'nAverages': n_avg,
                     'Tmax': tf,
                     'dt':   dt,
                     'pi2': pars['pi_half_len'],
                     'A_d':     amp_q_scaling,
                     'w_d':  pars['qb_LO']+qb_IF-detuning,
                     'w_LO': pars['qb_LO'],
            'wait_period':  resettime,
            }

    # save data
    with open(f"D:\weak_measurements\{exp}\data_{iteration:03d}.csv","w") as datafile:
        writer = csv.writer(datafile)
        writer.writerow(exp_dict.keys())
        writer.writerow(exp_dict.values())
        writer.writerow(t_arr)
        writer.writerow(I)
        writer.writerow(Q)

    return t_arr, I, Q,job,fitted_pars

def powerrabi(a_min = 0.01,    # minimum amp_q scaling
               a_max = 0.5,     # maximum amp_q scaling
               da = 0.005,       # step of amp_q
               n_avg = 2000,    # number of averages
               pulse_len = 1500,  # pulse length in nanoseconds
               plot = True):
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

    I, Q, job = getIQ(config, power_rabi, showprogress = True, n_avg = n_avg)

    if plot:
        pf.plot_data(x_vector=amps, y_vector=abs(I+1j*Q),sequence='a-rabi')

    return amps, I, Q;

#%% single_shot
def single_shot(nIterations=100000,
                n_reps = 1000,
                liveplot = False,
                numSamples = 1000,
                resetTime=400e3):

    resettime_clk = round(resetTime / 4)

    with program() as prog:
        i = declare(int) # number of iterations (used for continuously plotting)
        n = declare(int) # number of single shot measurements
        I = declare(fixed)
        Q = declare(fixed)
        Iexc = declare(fixed)
        Qexc = declare(fixed)
        I_st = declare_stream()
        Q_st = declare_stream()
        I_st_exc = declare_stream()
        Q_st_exc = declare_stream()
        i_st = declare_stream()
        # N_st = declare_stream()

        with for_(i, 0, i < nIterations, i + 1):
            with for_(n, 0, n < n_reps, n + 1):

                # do nothing
                align("qubit", "rr")
                measure("readout", "rr", None, *res_demod(I, Q))
                align('qubit','rr')
                wait(resettime_clk, "qubit")
                align('qubit','rr')
                # apply pi-pulse
                play("pi", "qubit")
                # play("gauss"*amp(0.41/0.45), "qubit",duration=round(274*2/4))
                align("qubit", "rr")
                measure("readout", "rr", None, *res_demod(Iexc, Qexc))
                wait(resettime_clk, "qubit")
                save(I, I_st)
                save(Q, Q_st)
                save(Iexc, I_st_exc)
                save(Qexc, Q_st_exc)
                # save(n,N_st)
            save(i, i_st)

        with stream_processing():
            # I_st.save('I')
            # Q_st.save('Q')
            # I_st_exc.save('I_exc')
            # Q_st_exc.save('Q_exc')
            # i_st.save('i')
            I_st.save_all('I')
            Q_st.save_all('Q')
            I_st_exc.save_all('Iexc')
            Q_st_exc.save_all('Qexc')
            i_st.save_all('i')

    datadict, job = get_results(config, prog,result_names=['I','Q','Iexc','Qexc', 'i'], showprogress=True, liveplot = liveplot)
    sourceFile = open('debug.py', 'w')
    print(generate_qua_script(debug_prog, config), file = sourceFile)
    sourceFile.close()
    return datadict, job,prog



#%% make_progress_meter
# display progress bar and send slack notification
def make_progress_meter(n_handle, n_total):

    # initialize counter
    n0 = 0

    while(n_handle.is_processing()):
        n_handle.wait_for_values(1)
    # create progressbar using tqdm
        with tqdm(n_handle.fetch_all(),total = n_total,miniters=1000,mininterval=1000) as progress_bar:
            pass
            # n = n_handle.fetch_all() # retrieve iteration value n
            # dn = n - n0 # calculate change in n since last update

            # if dn > 0:
            #     progress_bar.update(dn) # update progressbar with increase in n
            #     n0 = n # reset counter


# def plot_IQ_generic(I, Q, n, xdata = [0.0],  title = "", xlabel = ""):

#     fig, (ax1, ax2) = plt.subplots(2, sharex=True)
#     fig.suptitle(title, fontsize=11)

#     ax1.plot(xdata, mag(I, Q))
#     ax1.set_ylabel("mag (V)")

#     ax2.plot(xdata, phase(I, Q))
#     ax2.set_xlabel(xlabel)
#     ax2.set_ylabel("phase (rad)")

#     plt.tight_layout()

#     return fig;


#%% get_results
def get_results(config, jobtype, result_names = ["I", "Q", "N"],
                                showprogress = False,
                                nPoints = 1000,
                                notify = False,
                                liveplot = False):


    # Open Communication with the Server
    qmm = QuantumMachinesManager()

    # execute the job and get result handles
    qm = qmm.open_qm(config)
    job = qm.execute(jobtype)
    res_handles = job.result_handles


    # create a dictionary with result_names as keys
    # and their corresponding handles as values
    handles_dict = {}
    for name in result_names:
        handles_dict[name] = res_handles.get(name)

    # make live plot, if desired
    if liveplot:
        plot, ax = pf.init_IQ_plot()
        # wait for first values to come in before continuing to run python code
        for handle in handles_dict.values():
            handle.wait_for_values(1)
            is_processing = lambda: handle.is_processing()
        i0 = 0
        while(is_processing()):
            i = handles_dict["i"].fetch_all()['value'][-1] # retrieve iteration value n
            Δi = i - i0 # calculate change in n since last update
            if Δi > 0:
                datadict = get_data_from_handles(handles_dict)
                pf.plot_single_shot(datadict,axes=ax)
                # plot_IQ_blobs_plt(datadict, i)

    if showprogress:
        make_progress_meter(handles_dict['N'],n_total=nPoints)

    res_handles.wait_for_all_values()


    # retrieve all values
    datadict = get_data_from_handles(handles_dict)
    # plot_IQ_blobs_plt(datadict)

    # close quantum machine
    qmm.close_all_quantum_machines()

    return datadict, job;

def unpack_data(datadict, names = ["I", "Q"]):

    unpacked = map(lambda s: datadict[s], names)
    return list(unpacked)


def fix_length(*args):

    l = min(map(len, args))
    newargs = map(lambda a: a[0:l], args)
    return list(newargs)

#%% get_data_from_handles
def get_data_from_handles(handles_dict,verbose=0):

    datadict = {}

    for (key, handle) in handles_dict.items():
        # datadict[key] = handle.fetch_all()['value']
        datadict[key] = handle.fetch_all()
        if handle.has_dataloss():
            print(f'Dataloss occured in'+datadict[key])
        elif not handle.has_dataloss and verbose == 1:
            print('No data loss')

    return datadict


# datadict, job,debug_prog = single_shot(nIterations=10,resetTime=600e3,n_reps=1000, liveplot=True)


# #%% get IQ
# def getIQ(config, jobtype,  exp='rabi',    showprogress = True,
#                                 n_avg = 1000,
#                                 notify = False,
#                                 liveplot = False,
#                                 numSamples = 1000,
#                                 **plotkwargs):
#     """
#     getIQ(
#           config,       :: configuration dictionary (usually imported from a file named config.py),
#                            which is interpretted by QUA to determine settings on OPX
#           jobtype):     :: QUA jobtype (e.g. qubit_spec, rr_spec)
#         exp:            type of experiment running.
#     """

#     # Open Communication with the Server
#     qmm = QuantumMachinesManager()

#     # execute the job
#     qm = qmm.open_qm(config)
#     job = qm.execute(jobtype)
#     res_handles = job.result_handles

#     # get handles for I and Q
#     I_handle = res_handles.get("I")
#     Q_handle = res_handles.get("Q")

#     if exp == 'single_shot':
#         Iexc_handle = res_handles.get("Iexc")
#         Qexc_handle = res_handles.get("Qexc")

#     n_handle = res_handles.get('n')

#     # Wait until we know at least 1 value has been processed for this result (doesn't block subsequent python execution)
#     I_handle.wait_for_values(1)
#     Q_handle.wait_for_values(1)

#     if exp == 'single_shot':
#         Iexc_handle.wait_for_values(1)
#         Qexc_handle.wait_for_values(1)

#     n_handle.wait_for_values(1)

#     # display progress bar and send slack notification
#     if showprogress:

#         # get handle for iteration number n, and initialize counter
#         n0 = 0

#         # create progressbar using tqdm
#         with tqdm(total = n_avg) as progress_bar:

#             while(I_handle.is_processing()):

#                 # retrieve iteration value n
#                 n = n_handle.fetch_all()

#                 # send slack notification after 100 iterations
#                 if notify and n0 <= 100 and n > 100:

#                     later = gettimedelta(progress_bar)
#                     sendslack(message = f":zap: Measurement running. The measurement will be complete at approximately " + later + ".")

#                 # update progressbar with increase in n, and reset counter
#                 Δn = n - n0

#                 if Δn > 0:
#                     progress_bar.update(Δn)
#                     n0 = n

#     if liveplot:
#         assert exp == 'single_shot',"Live plotting only available for single shot!"
#         # make plot (initialize axes)
#         ax = pf.init_IQ_plot()


#         while(I_handle.is_processing()):
#             # wait until N measurements are completed
#             I_handle.wait_for_values(numSamples)
#             Q_handle.wait_for_values(numSamples)
#             Iexc_handle.wait_for_values(numSamples)
#             Qexc_handle.wait_for_values(numSamples)
#             # retrieve data
#             I = I_handle.fetch_all()
#             Q = Q_handle.fetch_all()
#             Iexc = Iexc_handle.fetch_all()
#             Qexc = Qexc_handle.fetch_all()

#             pf.plot_single_shot(I,Iexc,Q,Qexc,axes=ax)
#             job.resume()
#             # plotfunc(I, Q, n, **plotkwargs)
#     else:
#         res_handles.wait_for_all_values()

#     # retrieve I and Q values
#     I = I_handle.fetch_all()
#     Q = Q_handle.fetch_all()

#     # close quantum machine
#     qmm.close_all_quantum_machines()

#     return I, Q, job;

# def plot_IQ_blobs_init():

#     plot = sns.jointplot()
#     plot.set_axis_labels('I (mV)', 'Q (mV)')
#     plot.ax_marg_x.grid('on')
#     plot.ax_marg_y.grid('on')
#     plot.fig.tight_layout()
#     ax = plt.gca()
#     return plot, ax


# def plot_IQ_blobs_plt(datadict, i = -1):

#     I, Q, I_exc, Q_exc = fix_length(*unpack_data(datadict, names = ["I", "Q", "I_exc", "Q_exc"]))

#     ax = plt.gca()
#     ax.cla()

#     ax.scatter(I, Q, marker = 'o', c='orange', alpha=0.1)
#     ax.scatter(I_exc, Q_exc, marker = 'o', c='blue', alpha = 0.1)
#     if i >= 0:
#         plt.title(f"i = {i}")
#     plt.show()
#     plt.pause(0.1)
# def plot_IQ_blobs(plot, ax, datadict):

#     I = datadict["I"]
#     Q = datadict["Q"]
#     l = min(len(I), len(Q))
#     I = I[0:(l-1)]
#     Q = Q[0:(l-1)]

#     I_exc = datadict["I_exc"]
#     Q_exc = datadict["Q_exc"]
#     l_exc = min(len(I_exc), len(Q_exc))
#     I_exc = I[0:(l_exc-1)]
#     Q_exc = Q[0:(l_exc-1)]

#     plt.plot(I, Q)
#     plt.plot(I_exc, Q_exc)
#     plt.clf()

#     plot.ax_joint.plot(I*1e3, Q*1e3, 'o', label = "ground")
#     plot.ax_joint.plot(I_exc*1e3, Q_exc*1e3, 'o', label = "excited")
#     # plot.plot_marginals(sns.kdeplot)
#     # ax.legend()
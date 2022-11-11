"""
Created on Mon Oct 4 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
"""
from scipy.signal.windows import gaussian
from tqdm import tqdm
from qm import generate_qua_script
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
import plot_functions as pf
import os
from datetime import datetime
import csv
import glob
import time
import instrument_init as inst
import numpy as np
import json
from VISAdrivers.sa_api import *

class qubit():
#%% initialization
    def __init__(self,qb):
        try:
            with open(f'{qb}_pars.json', 'r') as openfile:
                self.pars = json.load(openfile)
            print('Loading parameter JSON file')
        except:
            print('Parameter file not found; loading parameters from template')
            self.pars = {
                "name":                         qb,
                "qubit_LO":                     int(3.0e9),
                "qubit_freq":                   int(4e9),
                "qubit_IF":                     50e6,
                "qubit_mixer_offsets":          [0,0], # I,Q
                "qubit_mixer_imbalance":        [0,0], # gain,phase
                "pi_half_len":                  260, # needs to be multiple of 4
                "pi_amp":                       0.44,
                "amp_q":                        0.45,
                "gauss_len":                    48,
                "gauss_amp":                    0.45,
                "rr_LO":                        int(6.55e9),
                "rr_freq":                      int(7),
                'rr_IF':                        50e6,
                "rr_mixer_offsets":             [0,0],
                "rr_mixer_imbalance":           [0,0],
                "amp_r":                        0.45,
                "tof":                          264, # time of flight in ns
                "rr_pulse_len_in_clk":          2500, # length of readout integration weights in clock cycles
                "IQ_rotation":                  -0/180*np.pi, # phase rotation applied to IQ data
                "analog_input_offsets":         [0,0]
                }

            with open(f'{qb}_pars.json', "w") as outfile:
                json.dump(self.pars, outfile)

        self.make_config(self.pars)

    #%% play_pulses
    def play_pulses(self):
        with program() as play_pulses:
            with infinite_loop_():
                play("const", 'qubit',duration=100)
                play("readout", "rr", duration=100)

        qmm = QuantumMachinesManager()
        qm = qmm.open_qm(self.config)
        job = qm.execute(play_pulses)

        return qm

    #%% tof_cal
    def tof_cal(self):
        qmm = QuantumMachinesManager()

        with program() as tof_cal:
            n = declare(int)
            adc_st = declare_stream(adc_trace=True)
            update_frequency('rr',10e6)
            with for_(n, 0, n < 1000, n + 1):
                reset_phase("rr")
                measure("readout", "rr", adc_st)
                wait(50000, "rr")
            with stream_processing():
                adc_st.input1().average().save("adc1")
                adc_st.input2().average().save("adc2")

        qm = qmm.open_qm(self.config)
        job = qm.execute(tof_cal)
        res_handles = job.result_handles
        res_handles.wait_for_all_values()
        adc1 = res_handles.get("adc1").fetch_all()
        adc2 = res_handles.get("adc2").fetch_all()

        pf.tof_plot(adc1, adc2)
        offset1 = np.mean(adc1)/4096
        offset2 = np.mean(adc2)/4096
        print(f'Input 1 Offset: {offset1*1e3} mV')
        print(f'Input 2 Offset: {offset2*1e3} mV')
        self.update_value('analog_input_offsets', value = [self.pars['analog_input_offsets'][0] - offset1,self.pars['analog_input_offsets'][1] - offset2])

    #%% punchout
    def punchout(self,
                 df = 0.1e6,
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
                dataI,dataQ,freqs,job = self.resonator_spec(f_LO=f_LO,atten=a,IF_min=df,IF_max=span,df=df,n_avg=n_avg,res_ringdown_time=res_ringdown_time,savedata=False,fit=False)
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
    def run_scan(self,
                 df = 0.1e6,
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
            numchunks = int((lo_max-lo_min)/chunksize) + 1
            lo_list = [i*chunksize+lo_min for i in range(numchunks)]
        else:
            numchunks = 1
            lo_list = [lo_min]

        for f in lo_list:
            if element == 'resonator':
                dataI,dataQ,freqs,job = self.resonator_spec(f_LO=f,atten=attenuation,IF_min=df,IF_max=chunksize,df=df,n_avg=n_avg,res_ringdown_time=res_ringdown_time,savedata=False)
            elif element == 'qubit':
                dataI,dataQ,freqs,job = self.qubit_spec(f_LO=f,amp_q_scaling=amp_q_scaling,saturation_dur=saturation_dur,atten=attenuation,IF_min=df,IF_max=chunksize,df=df,n_avg=n_avg,res_ringdown_time=res_ringdown_time,showprogress=showprogress,savedata=False)
            freq_arr.extend(freqs)
            I.extend(dataI)
            Q.extend(dataQ)

        exp_dict = {'date/time':    datetime.now(),
                   'nAverages': n_avg,
                         'w_LO': self.pars['rr_LO'],
                         'attenuation': attenuation,
                'wait_period':  res_ringdown_time,
                }

        # save data
        with open(f"D:\weak_measurements\spectroscopy\{element}_spec\data_{iteration:03d}.csv","w") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(exp_dict.keys())
            writer.writerow(exp_dict.values())
            writer.writerow(freq_arr)
            writer.writerow(I)
            writer.writerow(Q)

        return I, Q, freq_arr, job

    #%% resonator_spec
    def resonator_spec(self, IF_min = 0.1e6,
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

        self.update_value('rr_LO', value = f_LO)
        inst.set_rr_LO(self.pars['f_LO'])

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
                    measure("readout", "rr", None,*self.res_demod(I, Q))
                    save(I, I_st)
                    save(Q, Q_st)

                save(n,n_stream)

            with stream_processing():
                I_st.buffer(len(freqs)).average().save('I')
                Q_st.buffer(len(freqs)).average().save('Q')
                n_stream.save('n')

        datadict,job = self.get_results(rr_spec,result_names=["I","Q","n"],showprogress=False)

        I = datadict["I"]
        Q = datadict["Q"]
        freq_arr = freqs + self.pars['rr_LO']

        if fit:
            fc,fwhm = pf.fit_res(freq_arr,np.abs(I+1j*Q))
            pf.spec_plot(freq_arr,I,Q,attenuation=atten,df=df,iteration=iteration,element='resonator',fwhm=fwhm,fc=fc)
            print(f'Resonant Frequency: {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-6} MHz\nkappa = {2*np.pi*fwhm*1e-6:.3f} MHz')

        exp_dict = {'date/time':    datetime.now(),
                   'nAverages': n_avg,
                         'w_LO': self.pars['rr_LO'],
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
    def qubit_spec(self,
                   sa = 0,
                   f_LO = 5e9,
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

        try:
            list_of_files = glob.glob(r'D:\weak_measurements\spectroscopy\qubit_spec\*.csv')
            latest_file = max(list_of_files, key=os.path.getctime)
            iteration = int(latest_file[-7:-4].lstrip('0')) + 1
        except:
            iteration = 1

        freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
        # freqs_list = freqs.tolist()
        saturation_dur = int(saturation_dur)
        inst.set_attenuator(attenuation=atten)


        self.update_value('qubit_LO',value = f_LO)
        self.update_value('rr_LO', value = self.pars['rr_freq'] - self.pars['rr_IF'])

        self.check_mix_cal(sa,check=True,threshold = - 50)
        # self.opt_mixer(sa,cal='SB',mode='coarse',element='qubit',reference=-30)
        # self.opt_mixer(sa,cal='SB',mode='fine',element='qubit',reference=-50)

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
                # with for_each_(f, freqs_list):
                    # update IF frequency going into qubit mixer
                    update_frequency("qubit", f)
                    # measure background
                    if on_off:
                        measure("readout", "rr", None, *self.res_demod(I_background, Q_background))
                        wait(res_ringdown_time, "rr")
                        align("rr", "qubit") # wait for operations on resonator to finish before playing qubit pulse
                    # play qubit pulse and measure
                    play("const" * amp(amp_q_scaling), "qubit", duration = saturation_dur)
                    align("qubit", "rr") # wait for operations on resonator to finish before playing qubit pulse
                    measure("readout", "rr", None, *self.res_demod(I, Q))
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
        datadict, job = self.get_results(QubitSpecProg, result_names = ["I", "Q", "N"], nPoints=n_avg,showprogress=showprogress, notify = notify)
        I = datadict["I"]
        Q = datadict["Q"]
        freq_arr = freqs + self.pars['qubit_LO']

        pf.spec_plot(freq_arr,I,Q,iteration=iteration,element='qubit',find_peaks=True)
        # print(f'Qubit Frequency: {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-6} MHz\nkappa = {2*np.pi*fwhm*1e-6:.3f} MHz')

        if savedata:
            exp_dict = {'date/time':     datetime.now(),
                       'nAverages': n_avg,
                             'A_d':     amp_q_scaling,
                             'w_LO': self.pars['qubit_LO'],
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

        return I, Q, freq_arr, job;

    #%% pulse_exp
    def pulse_exp(self,exp='rabi',
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
                    with for_(t, t_min, t < tmax + dt/2, t + dt):  # Sweep pulse duration

                        play("gauss" * amp(amp_q_scaling), "qubit", duration=t)
                        align("qubit", "rr")
                        measure("readout", "rr", None, *self.res_demod(I, Q))
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
                           measure("readout", "rr", None,*self.res_mod(I,Q))
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
                           measure("readout", "rr", None, *self.res_demod(I, Q))
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
                            measure("readout", "rr", None,*self.res_mod(I,Q))
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
                            measure("readout", "rr", None, *self.res_demod(I, Q))
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
                            measure("readout", "rr", None, *self.res_mod(I,Q))
                            save(I, I_stream)
                            save(Q, Q_stream)

                        with else_():

                            wait(resettime_clk, "qubit")
                            play("pi", "qubit")
                            align("qubit", "rr")
                            wait(t, 'rr')
                            measure("readout", "rr", None,*self.res_mod(I,Q))
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
                         'w_d':  pars['qubit_LO']+qb_IF-detuning,
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

    #%% get_results
    def get_results(self,jobtype, result_names = ["I", "Q", "n"],
                                    showprogress = False,
                                    nPoints = 1000,
                                    notify = False,
                                    liveplot = False):
        """

        Args:
            jobtype (QUA program): QUA program to execute.
            result_names (array of strings, optional): The names of the variables we want to stream. Defaults to ["I", "Q", "n"].
            showprogress (TYPE, optional): DESCRIPTION. Defaults to False.
            nPoints (TYPE, optional): DESCRIPTION. Defaults to 1000.
            notify (TYPE, optional): DESCRIPTION. Defaults to False.
            liveplot (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            datadict (TYPE): DESCRIPTION.
            TYPE: DESCRIPTION.

        """


        # Open Communication with the Server
        qmm = QuantumMachinesManager()

        # execute the job and get result handles
        qm = qmm.open_qm(self.config)
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
                    datadict = self.get_data_from_handles(handles_dict)
                    pf.plot_single_shot(datadict,axes=ax)
                    # plot_IQ_blobs_plt(datadict, i)

        if showprogress:
            make_progress_meter(handles_dict['N'],n_total=nPoints)

        res_handles.wait_for_all_values()


        # retrieve all values
        datadict = self.get_data_from_handles(handles_dict)
        # plot_IQ_blobs_plt(datadict)

        # close quantum machine
        qmm.close_all_quantum_machines()

        return datadict, job;

    def unpack_data(self,datadict, names = ["I", "Q"]):

        unpacked = map(lambda s: datadict[s], names)
        return list(unpacked)


    def fix_length(*args):

        l = min(map(len, args))
        newargs = map(lambda a: a[0:l], args)
        return list(newargs)

    #%% get_data_from_handles
    def get_data_from_handles(self,handles_dict,verbose=0):

        datadict = {}

        for (key, handle) in handles_dict.items():
            # datadict[key] = handle.fetch_all()['value']
            datadict[key] = handle.fetch_all()
            if handle.has_dataloss():
                print(f'Dataloss occured in'+datadict[key])
            elif not handle.has_dataloss and verbose == 1:
                print('No data loss')

        return datadict

    #%% config_sa
    def config_sa(self,sa,freq,span=5e6,reference=-30):
        """
        Prepares spectrum analyzer for measurement

        Parameters
        ----------
        sa :
            Handle for spectrum analyzer
        freq : float
            Center frequency of span.
        span : float, optional
            DESCRIPTION. The default is 5e6.
        reference : float, optional
            Upper power threshold of SA in dBm. The default is -30.

        Returns
        -------
        None.

        """

        sa_config_level(sa, reference) # sets sensitivity
        sa_config_center_span(sa, freq, span) # sets center frequency
        sa_initiate(sa, SA_SWEEPING, 0)
        query = sa_query_sweep_info(sa)
        sweep_length = query["sweep_length"]
        start_freq = query["start_freq"]
        bin_size = query["bin_size"]
        freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)

    #%% get_power
    def get_power(self,sa,freq,reference=-100,span=1e6, config=False, plot=False, output = True):
        """
        Configures SA (optional) and measures power at specified frequency

        Parameters
        ----------
        sa : ???
            spectrum analyzer.
        freq : TYPE
            DESCRIPTION.
        reference : TYPE, optional
            DESCRIPTION. The default is -100.
        span : TYPE, optional
            DESCRIPTION. The default is 5e6.
        config : boolean, optional
            whether to reconfigure the SA or not. Set to false when calibrating mixers. The default is False.
        plot : TYPE, optional
            DESCRIPTION. The default is False.
        output: boolean, optional
            whether or not to print the power at the requested frequency

        Returns
        -------
        freqs : TYPE
            DESCRIPTION.
        power : TYPE
            DESCRIPTION.

        """
        atten = inst.get_attenuation()
        inst.set_attenuator(attenuation=0)

        # skips configuring the spectrum analyzer. Used only when optimizing mixer
        if config:
            self.play_pulses()
            sa_config_level(sa, reference) # sets sensitivity
            sa_config_center_span(sa, freq, span) # sets center frequency
            sa_initiate(sa, SA_SWEEPING, 0)
            query = sa_query_sweep_info(sa)
            sweep_length = query["sweep_length"]
            start_freq = query["start_freq"]
            bin_size = query["bin_size"]
            freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)

        # measure
        signal = sa_get_sweep_64f(sa)['max']
        power = round(np.max(signal),1)

        if plot:
            pf.power_plot(freqs, signal, power, fc=freq)
        if output:
            print(f'{power} dBm at {freq/1e9} GHz')
        inst.set_attenuator(attenuation=atten)
        return power

    #%% opt_mixer
    def opt_mixer(self,sa,cal,mode,element,freq_span=1e6,reference = -30, plot=True):
        """
        Minimizes leakage at LO ('lo' option) or at image sideband ('sb' option) by sweeping the relevant parameters

        Args:
            sa (): spectrum analyzer handle.
            cal (str): Whether to minimize LO leakage or image sideband. The default is 'lo'.
            mode (str): Coarse of fine stepsize.
            element (str):       Which element to optimize (qubit (qubit) or readout (rr)).
            pars (dict):        dictionary containing experimental values like mixer calibration offsets.
            reference (float): Threshold of spectrum analyzer.
            plot (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            values (float array): optimal calibration values.
            argmin (TYPE): DESCRIPTION.

        """
        qm = self.play_pulses()
        # gets frequency values from config file
        freqLO = self.config['mixers'][f'{element}'][0]['lo_frequency']
        freqIF = self.config['mixers'][f'{element}'][0]['intermediate_frequency']

        if cal == 'LO':
            freq = freqLO
            par1 = self.pars[f'{element}_mixer_offsets'][0]
            par2 = self.pars[f'{element}_mixer_offsets'][1]
            print(f'LO at {round(freqLO*1e-9,5)} GHz\nCurrent I_offset = {round(par1*1e3,1)} mV, Current Q_offset = {round(par2*1e3,1)} mV')
        elif cal == 'SB':
            freq = freqLO - freqIF
            par1 = self.pars[f'{element}_mixer_imbalance'][0]
            par2 = self.pars[f'{element}_mixer_imbalance'][1]
            # if np.abs(par2) > 150e-3:
            #     par2 = 0
            print(f'Sideband at {round((freqLO-freqIF)*1e-9,5)} GHz\nCurrent gain = {round(par1,4)}, Current phase = {round(par2,4)}')

        if element == 'rr':
            atten = inst.get_attenuation()
            inst.set_attenuator(attenuation=0)

        self.config_sa(sa,freq,span=freq_span,reference=reference) # setup spectrum analyzer for measurement

        # initialize sweep parameters
        if cal == 'LO':
            if mode == 'coarse':
                span = 40e-3
                step = 10e-3
            elif mode == 'intermediate':
                span = 10e-3
                step = 2e-3
            elif mode == 'fine':
                span = 1e-3
                step = 0.1e-3
        elif cal == 'SB':
            if mode == 'coarse':
                span = 150e-3
                step = 30e-3
            elif mode == 'intermediate':
                span = 50e-3
                step = 10e-3
            elif mode == 'fine':
                span = 10e-3
                step = 1e-3

        par1_arr = np.arange(par1-span/2, par1+span/2, step)
        par2_arr = np.arange(par2-span/2, par2+span/2, step)
        L1 = len(par1_arr)
        L2 = len(par2_arr)
        power_data = np.zeros((L1,L2))

        # sweep parameters and get power at every point
        with tqdm(total = L1*L2) as progress_bar:
            for i, par1 in enumerate((par1_arr)):
                for j, par2 in enumerate((par2_arr)):
                    if cal == 'LO':
                        qm.set_output_dc_offset_by_element(element, "I", par1)
                        qm.set_output_dc_offset_by_element(element, "Q", par2)
                    elif cal == 'SB':
                        qm.set_mixer_correction(element,int(freqIF), int(freqLO), self.IQ_imbalance(par1, par2))
                    time.sleep(0.1)
                    power_data[i,j] = self.get_power(sa, freq,span=freq_span,output=False)
                    progress_bar.update(1)

        argmin = np.unravel_index(np.argmin(power_data), power_data.shape)
        # print(argmin)
        # set the parameters to the optimal values and modify the JSON dictionary
        if cal == 'LO':
            opt_I = par1_arr[argmin[0]]
            opt_Q = par2_arr[argmin[1]]
            qm.set_output_dc_offset_by_element(element, "I", opt_I)
            qm.set_output_dc_offset_by_element(element, "Q", opt_Q)
            self.update_value(f'{element}_mixer_offsets',[opt_I,opt_Q])
            print(f'optimal I_offset = {round(opt_I*1e3,1)} mV, optimal Q_offset = {round(1e3*opt_Q,1)} mV')
        elif cal == 'SB':
            qm.set_mixer_correction(element,int(freqIF), int(freqLO), self.IQ_imbalance(par1_arr[argmin[0]],par2_arr[argmin[1]]))
            self.update_value(f'{element}_mixer_imbalance',[par1_arr[argmin[0]],par2_arr[argmin[1]]])
            print(f"optimal gain = {round(par1_arr[argmin[0]],4)}, optimal phi = {round(par2_arr[argmin[1]],4)}")

        if element == 'rr':
            inst.set_attenuator(atten)

        print(f'Power: {np.amin(power_data)} dBm at {freq/1e9} GHz')
        if plot:
            pf.plot_mixer_opt(par1_arr, par2_arr, power_data,cal=cal,element=element,fc=freq)

    #%% check_mixer_calibration
    def check_mix_cal(self, sa, check = True, threshold = -50):
        # checks LO leakage and image sideband power and opt mixer accordingly
        if check:
            for element in ['qubit','rr']:
                for j in ['LO','SB']:
                    if j == 'LO':
                        fc = self.pars[f'{element}_LO']
                    elif j == 'SB':
                        fc = self.pars[f'{element}_LO'] -  self.pars[f'{element}_IF']
                    print(f'Checking {element} {j}')
                    leak = self.get_power(sa, freq = fc, reference = -20, config = True, plot = False)
                    while leak > threshold:
                        leak0 = leak
                        ref = leak + 10
                        print(f'Minimizing {element} {j} leakage')
                        if leak > - 40:
                            self.opt_mixer(sa,cal = j, mode = 'coarse', element = element, reference = ref)
                        elif leak < - 40 and leak > - 55:
                            self.opt_mixer(sa, cal = j, mode = 'intermediate', element = element, reference = ref)
                        elif leak < - 55:
                            self.opt_mixer(sa, cal = j, mode = 'fine', element = element, reference = ref)

                        leak = self.get_power(sa,freq = fc,reference = ref, config = True, plot = False)

                        if (leak - leak0) < 3:
                            print("Can't optimize mixer further")
                            break
        else:
            pass

    #%% meas_utils
    def res_demod(self,I, Q):

        return (dual_demod.full("integW_cos", "out1", "integW_minus_sin", "out2", I),
                dual_demod.full("integW_sin", "out1", "integW_cos", "out2", Q))

    def IQ_imbalance(self,g, phi):
        c = np.cos(phi)
        s = np.sin(phi)
        N = 1 / ((1-g**2)*(2*c**2-1))
        return [float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]

    def delayed_gauss(self,amp=0.45, length=4, sigma=2):
        gauss_arg = np.linspace(-sigma, sigma, length)
        delay = 16 - length - 4
        if delay < 0:
            return amp * np.exp(-(gauss_arg ** 2) / 2)

        return np.r_[np.zeros(delay), amp * np.exp(-(gauss_arg ** 2) / 2), np.zeros(4)]

    def update_value(self,key,value):
        self.pars[key] = value
        self.make_config(self.pars)

        if key == 'qubit_LO':
            inst.set_qb_LO(value)
        elif key == 'rr_LO':
            inst.set_rr_LO(value)

        with open(f'{self.pars["name"]}_pars.json', "w") as outfile:
            json.dump(self.pars, outfile)

    #%% make_configuration_file
    def make_config(self,pars):

        gauss_wf_4ns = self.delayed_gauss()

        self.config = {

            "version": 1,

            "controllers": {
                "con1": {
                    "type": "opx1",
                    "analog_outputs": {
                        1: {"offset": pars['qubit_mixer_offsets'][0]},  # qubit I
                        2: {"offset": pars['qubit_mixer_offsets'][1]},  # qubit Q
                        3: {"offset": pars['rr_mixer_offsets'][0]},  # rr I
                        4: {"offset": pars['rr_mixer_offsets'][1]},  # rr Q
                    },
                    "digital_outputs": {},
                    "analog_inputs": {
                        1: {"offset": pars['analog_input_offsets'][0], "gain_db": 3},  # rr I
                        2: {"offset": pars['analog_input_offsets'][1], "gain_db": 3}  # rr Q
                    },
                },
            },

            "elements": {
                "qubit": {
                    "mixInputs": {
                        "I": ("con1", 1),
                        "Q": ("con1", 2),
                        "lo_frequency": pars['qubit_LO'],
                        "mixer": "qubit",
                    },
                    "intermediate_frequency": 50e6,
                    "digitalInputs": {},
                    "operations": {
                        "const": "const_pulse_IQ",
                        "gauss": "gaussian_pulse",
                        "gauss_4ns": "gaussian_4ns_pulse",
                        "pi": "pi_pulse1",
                        "pi_half": "pi_half_pulse1",
                        "arb_op": "arb_pulse",
                        "X": "Xpi_pulse",
                        "Y": "Ypi_pulse",
                        "X/2": "Xpi_half_pulse",
                        "Y/2": "Ypi_half_pulse",
                    },
                },
                "rr": {
                    "mixInputs": {
                        "I": ("con1", 3),
                        "Q": ("con1", 4),
                        "lo_frequency": pars['rr_LO'],
                        "mixer": "rr",
                    },
                    "intermediate_frequency": 50e6,
                    "outputs": {
                        "out1": ("con1", 1),
                        "out2": ("con1", 2),
                    },
                    "time_of_flight": pars['tof'], # should be multiple of 4 (at least 24)
                    "smearing": 0, # adds 40ns of data from each side of raw adc trace to account for ramp up and down of readout pulse
                    "operations": {
                        "const": "const_pulse_IQ_rr",
                        "readout": "ro_pulse1",
                    },
                },
            },

            "pulses": {
                "const_pulse_IQ": {
                    "operation": "control",
                    "length": 100,
                    "waveforms": {
                        "I": "const_wf",
                        "Q": "zero_wf",
                    },
                },
                "const_pulse_IQ_rr": {
                    "operation": "control",
                    "length": 100,
                    "waveforms": {
                        "I": "const_wf_rr",
                        "Q": "zero_wf",
                    },
                },
                "pi_pulse1": {
                    "operation": "control",
                    "length": 2*pars['pi_half_len'],
                    "waveforms": {
                        "I": "pi_wf_i1",
                        "Q": "pi_wf_q1",
                    },
                },
                "Xpi_pulse": {
                    "operation": "control",
                    "length": 2*pars['pi_half_len'],
                    "waveforms": {
                        "I": "pi_wf_i1",
                        "Q": "pi_wf_q1",
                    },
                },
                "Ypi_pulse": {
                    "operation": "control",
                    "length": 2*pars['pi_half_len'],
                    "waveforms": {
                        "I": "pi_wf_q1",
                        "Q": "pi_wf_i1",
                    },
                },
                "Xpi_half_pulse": {
                    "operation": "control",
                    "length": pars['pi_half_len'],
                    "waveforms": {
                        "I": "pi_half_wf_i1",
                        "Q": "pi_half_wf_q1",
                    },
                },
                "Ypi_half_pulse": {
                    "operation": "control",
                    "length": pars['pi_half_len'],
                    "waveforms": {
                        "I": "pi_half_wf_q1",
                        "Q": "pi_half_wf_i1",
                    },
                },
                "gaussian_pulse": {
                    "operation": "control",
                    "length": pars['gauss_len'],
                    "waveforms": {
                        "I": "gaussian_wf",
                        "Q": "zero_wf",
                    },
                },
                "gaussian_4ns_pulse": {
                    "operation": "control",
                    "length": 16,
                    "waveforms": {
                        "I": "gaussian_4ns_wf",
                        "Q": "zero_wf",
                    },
                },
                "pi_half_pulse1": {
                    "operation": "control",
                    "length": pars['pi_half_len'],
                    "waveforms": {
                        "I": "pi_half_wf_i1",
                        "Q": "pi_half_wf_q1",
                    },
                },
                "ro_pulse1": {
                    "operation": "measurement",
                    "length": pars['rr_pulse_len_in_clk'] *4, # in ns (needs to be multiple of 4)
                    "waveforms": {"I": "ro_wf1", "Q": "zero_wf"},
                    "integration_weights": {
                        "integW_cos": "integW1_cos",
                        "integW_sin": "integW1_sin",
                        "integW_minus_sin": "integW1_minus_sin"
                    },
                    "digital_marker": "ON",
                },
                "arb_pulse": {
                    "operation": "control",
                    "length": 40,
                    "waveforms": {
                        "I": "arb_wfm",
                        "Q": "zero_wf",
                    },
                },
            },

            "waveforms": {
                "zero_wf": {"type": "constant", "sample": 0.0},
                "const_wf": {"type": "constant", "sample": pars['amp_q']},
                "const_wf_rr": {"type": "constant", "sample": pars['amp_r']},
                "gaussian_wf": {"type": "arbitrary", "samples": [float(arg) for arg in pars['gauss_amp'] * gaussian(pars['gauss_len'], pars['gauss_len']/5)]},
                "gaussian_4ns_wf": {"type": "arbitrary", "samples": gauss_wf_4ns},
                "ro_wf1": {"type": "constant", "sample": pars['amp_r']},
                "pi_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_amp'] * gaussian(2*pars['pi_half_len'], 2*pars['pi_half_len']/5)]},
                "pi_wf_q1": {"type": "constant", "sample": 0.0},
                "pi_half_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_amp'] * gaussian(pars['pi_half_len'], pars['pi_half_len']/5)]},
                "pi_half_wf_q1": {"type": "constant", "sample": 0.0},
                "arb_wfm": {"type": "arbitrary", "samples": [0.2]*10+[0.3]*10+[0.25]*20},
            },

            "digital_waveforms": {
                "ON": {"samples": [(1, 0)]}
            },

            "integration_weights": {
                "integW1_cos": {
                    "cosine": [(np.cos(pars['IQ_rotation']) , pars['rr_pulse_len_in_clk']*4)],
                    "sine": [(-np.sin(pars['IQ_rotation']) , pars['rr_pulse_len_in_clk']*4)],
                },
                "integW1_sin": {
                    "cosine": [(np.sin(pars['IQ_rotation']) , pars['rr_pulse_len_in_clk']*4)],
                    "sine": [(np.cos(pars['IQ_rotation']) , pars['rr_pulse_len_in_clk']*4)],
                },
                "integW1_minus_sin": {
                    "cosine": [(-np.sin(pars['IQ_rotation']) , pars['rr_pulse_len_in_clk']*4)],
                    "sine": [(-np.cos(pars['IQ_rotation']) , pars['rr_pulse_len_in_clk']*4)],
                },
                "integW2_cos": {
                    "cosine": [1.0] * pars['rr_pulse_len_in_clk'],
                    "sine": [0.0] * pars['rr_pulse_len_in_clk'],
                },
                "integW2_sin": {
                    "cosine": [0.0] * pars['rr_pulse_len_in_clk'],
                    "sine": [1.0] * pars['rr_pulse_len_in_clk'],
                },
                "integW2_minus_sin": {
                    "cosine": [0.0] * pars['rr_pulse_len_in_clk'],
                    "sine": [-1.0] * pars['rr_pulse_len_in_clk'],
                }
            },

            "mixers": {
                "qubit": [{"intermediate_frequency": 50e6, "lo_frequency": pars['qubit_LO'], "correction": self.IQ_imbalance(*pars['qubit_mixer_imbalance'])}],
                "rr": [{"intermediate_frequency": 50e6, "lo_frequency": pars['rr_LO'], "correction": self.IQ_imbalance(*pars['rr_mixer_imbalance'])}],
            }
        }

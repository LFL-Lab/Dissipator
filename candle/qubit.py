"""
Created on Mon Oct 4 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
"""
from scipy.signal.windows import gaussian
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
from qualang_tools.loops import from_array
from qualang_tools.analysis.discriminator import two_state_discriminator
from qm.logger import logger
from Utilities import *

host='10.71.0.56'
port='9510'
logger.setLevel(level='WARNING')

class qubit():
#%% INITIALIZATION

    #%%% default_pars
    default_pars = {
                    "qubit_LO":                     int(4.48e9),
                    "qubit_freq":                   int(4.5129e9),
                    "qubit_IF":                     int(4.5129e9) - int(4.48e9),
                    "qubit_mixer_offsets":          [0,0], # I,Q
                    "qubit_mixer_imbalance":        [0,0], # gain,phase
                    "pi_len":                       48, # needs to be multiple of 4
                    "pi_half_len":                  48, # needs to be multiple of 4
                    "pi_half_amp":                  0.2,
                    "pi_amp":                       0.45,
                    "amp_q":                        0.45,
                    "gauss_len":                    48,
                    "gauss_amp":                    0.45,
                    "rr_LO":                        int(6.42e9),
                    "rr_freq":                      int(6.4749e9),
                    'rr_IF':                        int(6.4749e9) - int(6.42e9),
                    "rr_mixer_offsets":             [0,0],
                    "rr_mixer_imbalance":           [0,0],
                    "amp_r":                        0.45,
                    'rr_atten':                     25,
                    "tof":                          260, # time of flight in ns
                    "rr_pulse_len_in_clk":          500, # length of readout integration weights in clock cycles
                    "IQ_rotation":                  -0/180*np.pi, # phase rotation applied to IQ data
                    "analog_input_offsets":         [0,0],
                    "qubit_resettime":              400e3,
                    "rr_resettime":                 20e3
                    }

    #%%% __init__
    def __init__(self, qb):
        # load pars from json, OR create new json file
        self.name = qb
        try:
            print('Loading parameter JSON file')
            with open(f'{qb}_pars.json', 'r') as openfile:
                self.pars = json.load(openfile)

            # compare keys
            default_keys = set(self.default_pars.keys())
            keys = set(self.pars.keys())

            # find all the keys in default_pars that are not in pars, and add them to pars w/ default value
            for k in (default_keys - keys):
                self.add_key(k, self.default_pars[k])

            # find all the keys in pars that are not in default_pars, and remove them from pars
            for k in (keys - default_keys):
                self.remove_key(k)


        except FileNotFoundError:
            print('Parameter file not found; loading parameters from template')
            self.pars = self.default_pars


        self.write_pars()
        self.init_instruments()
        self.make_config(self.pars)


#%% EXPERIMENTS
 #%%% play_pulses
    def play_pulses(self,amp_q=1):
        with program() as play_pulses:
            with infinite_loop_():
                play("const"*amp(amp_q), 'qubit',duration=100)
                play("readout", "rr", duration=100)

        qmm = QuantumMachinesManager(host=host, port=port)
        qm = qmm.open_qm(self.config)
        job = qm.execute(play_pulses)

        return qm


    #%%% punchout
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


    #%%% run_scan
    def run_scan(self,
                 df = 0.1e6,
                 n_avg = 500,
                 element='resonator',
                 check_mixers=True,
                 chunksize = 200e6,
                 attenuation=20,
                 lo_min = 6e9,
                 lo_max = 7e9,
                 amp_q_scaling = 1,
                 saturation_dur = 20e3,
                 showprogress=False):
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
                dataI,dataQ,freqs,job = self.resonator_spec(f_LO=f,atten=attenuation,check_mixers=check_mixers,IF_min=df,IF_max=chunksize,df=df,n_avg=n_avg,res_ringdown_time=res_ringdown_time,savedata=False)
            elif element == 'qubit':
                dataI,dataQ,freqs,job = self.qubit_spec(f_LO=f,amp_q_scaling=amp_q_scaling,check_mixers=check_mixers,saturation_dur=saturation_dur,atten=attenuation,IF_min=df,IF_max=chunksize,df=df,n_avg=n_avg,showprogress=showprogress,savedata=False)
            freq_arr.extend(freqs)
            I.extend(dataI)
            Q.extend(dataQ)

        exp_dict = {'date/time':    datetime.now(),
                   'nAverages': n_avg,
                         'w_LO': self.pars['rr_LO'],
                         'attenuation': attenuation,
                'wait_period':  self.pars['rr_ringdown_time'],
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

    #%%% resonator_spec
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

        I = np.array(datadict["I"])
        Q = np.array(datadict["Q"])
        freq_arr = np.array(freqs + self.pars['rr_LO'])

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

    #%%% qubit_spec
    def qubit_spec(self,
                   sa = 0,
                   f_LO = 5e9,
                   IF_min = 0.1e6,          # min IF frequency
                   IF_max = 400e6,          # max IF frequency
                   check_mixers=True,
                   df = 0.1e6,              # IF frequency step
                   rr_freq = 6e9,       #resonator frequency
                   amp_q_scaling = 0.1,     # prefactor to scale default "const" qubit tone, amp_q
                   n_avg = 500, # number of averages
                   atten = 10, # readout attenuation
                   saturation_dur = int(10e3),   # time qubit saturated w/ qubit tone, in ns
                   resettime = int(40e3),      # wait time between experiments, in ns
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

        freq_arr = np.arange(IF_min, IF_max + df/2, df, dtype=int)
        saturation_dur = int(saturation_dur)
        inst.set_attenuator(attenuation=atten)


        self.update_value('qubit_LO',value = f_LO)
        self.update_value('rr_LO', value = self.pars['rr_freq'] - self.pars['rr_IF'])

        self.check_mix_cal(sa,check=check_mixers,amp_q = amp_q_scaling, threshold = - 60)

        prog = self.make_sequence(exp='qubit-spec',on_off=on_off,saturation_dur=saturation_dur,var_arr=freq_arr,n_avg=n_avg,amp_q_scaling=amp_q_scaling,resettime_clk=clk(resettime))
        # execute 'QubitSpecProg' using configuration settings in 'config'
        # fetch averaged I and Q values that were saved
        datadict, job = self.get_results(jobtype = prog, result_names = ["I", "Q"], n_total=n_avg,showprogress=showprogress, notify = notify)

        qb_power = self.get_power(sa,freq=self.pars['qubit_LO']+self.pars['qubit_IF'],reference=0,amp_q = amp_q_scaling, span=1e6,config=True,output=False)
        rr_power = self.get_power(sa,freq=self.pars['rr_LO']+self.pars['rr_IF'],reference=0,span=1e6,config=True,output=False)

        I = np.array(datadict["I"])
        Q = np.array(datadict["Q"])
        freq_arr = np.array(freq_arr + self.pars['qubit_LO'])

        pf.spec_plot(freq_arr,I,Q,iteration=iteration,element='qubit',qb_power=qb_power,rr_power=rr_power,rrFreq=self.pars['rr_freq'],find_peaks=True)
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


    #%%% power_rabi
    def power_rabi(self,sa = 0,
                   a_min = 0.01,    # minimum amp_q scaling
                   a_max = 1,     # maximum amp_q scaling
                   da = 0.005,       # step of amp_q
                   check_mixers = True,
                   pulse = 'pi',
                   n_avg = 2000,    # number of averages
                   fit = True,
                   plot = True,
                   detuning = 0e6):

        amps = np.arange(a_min, a_max + da/2, da)

        inst.set_attenuator(attenuation=self.pars['rr_atten'])

        self.check_mix_cal(sa, check = check_mixers, threshold = -55)

        prog = self.make_sequence(exp = 'p-rabi', pulse = pulse,
                                  var_arr = amps,
                                  detuning = detuning,
                                  n_avg = n_avg)

        datadict, job = self.get_results(jobtype = prog, result_names = ['I','Q'], showprogress = True, progress_key = 'n', n_total = n_avg)

        I = datadict['I']
        Q = datadict['Q']

        if fit:
            fitted_pars,error = pf.fit_data(amps,np.abs(I+1j*Q),sequence='p-rabi',dt=amps[-1]/len(amps),fitFunc='rabi')
            if plot:
                pf.plot_data(x_vector=amps, y_vector=np.abs(I+1j*Q),sequence='p-rabi',fitted_pars=fitted_pars,
                             qubitDriveFreq=self.pars['qubit_LO']+self.pars['qubit_IF'],savefig=False,nAverages=n_avg)

        '''Update pulse amplitude'''
        A = fitted_pars[1] #self.pars['pi_amp'] * fitted_pars[1]

        return amps, I, Q, job, A


    #%%% single_shot
    def single_shot(self,nIterations=100000,
                    n_reps = 1000,
                    liveplot = False,
                    numSamples = 1000):

        inst.set_attenuator(attenuation=self.pars['rr_atten'])

        prog = self.make_sequence(exp='ss', nIterations = nIterations, n_reps = n_reps)

        datadict, job = self.get_results(jobtype = prog,result_names=['i', 'I','Q','Iexc','Qexc'], showprogress=False, liveplot = liveplot)

        if liveplot == False:
            plot, ax = pf.init_IQ_plot()
            pf.plot_single_shot(datadict,axes=ax)

        return datadict, job, prog

    #%%% make_sequence
    def make_sequence(self,exp='rabi',var_arr=0,detuning=0,n_avg=0,amp_q_scaling=1,numPeriods=2,nIterations=1,n_reps=100,on_off=True,saturation_dur=int(10e3),pulse='pi'):

        resettime_clk= clk(self.pars['qubit_resettime'])

        if exp == 'qubit-spec':
            with program() as prog:
                n, f, I, Q = self.declare_vars([int, int, fixed, fixed])

                I_stream, Q_stream, n_stream = self.declare_streams(stream_num=3)

                if on_off:
                    I_b, Q_b, I_tot,Q_tot = self.declare_vars([fixed,fixed, fixed, fixed])

                # loop over n_avg iterations
                with for_(n, 0, n < n_avg, n + 1):
                    save(n,n_stream)
                    # loop over list of IF frequencies
                    with for_(*from_array(f,var_arr)):
                        # update IF frequency going into qubit mixer
                        update_frequency("qubit", f)
                        # measure background
                        if on_off:
                            measure("readout", "rr", None, *self.res_demod(I_b, Q_b))
                            wait(clk(self.pars['rr_resettime']), "rr")
                            align("rr", "qubit") # wait for operations on resonator to finish before playing qubit pulse
                        # play qubit pulse and measure
                        play("const" * amp(amp_q_scaling), "qubit", duration = clk(saturation_dur))
                        align("qubit", "rr") # wait for operations on resonator to finish before playing qubit pulse
                        measure("readout", "rr", None, *self.res_demod(I, Q))
                        # subtract background and save to stream
                        if on_off:
                            assign(I_tot, I - I_b)
                            assign(Q_tot, Q - Q_b)
                            save(I_tot, I_stream)
                            save(Q_tot, Q_stream)
                        else:
                            save(I, I_stream)
                            save(Q, Q_stream)
                        # wait some time before continuing to next IF frequency
                        wait(resettime_clk, "rr")

                # average data over iterations and save to stream
                with stream_processing():
                    I_stream.buffer(len(var_arr)).average().save('I')
                    Q_stream.buffer(len(var_arr)).average().save('Q')
                    n_stream.save('n')

        if exp == 'rabi':
            with program() as prog:

                n, t, I, Q = self.declare_vars([int, int, fixed, fixed])

                I_stream, Q_stream, n_stream = self.declare_streams(stream_num=3)

                update_frequency('qubit', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning) # sets the IF frequency of the qubit

                with for_(n, 0, n < n_avg, n + 1):
                    save(n, n_stream)
                    with for_(*from_array(t,var_arr)):
                        play("gauss" * amp(amp_q_scaling), "qubit", duration=t)
                        align("qubit", "rr")
                        measure("readout", "rr", None, *self.res_demod(I, Q))
                        # save(t,t_stream)
                        save(I, I_stream)
                        save(Q, Q_stream)
                        wait(resettime_clk,"qubit")

                with stream_processing():
                    I_stream.buffer(len(var_arr)).average().save("I")
                    Q_stream.buffer(len(var_arr)).average().save("Q")
                    n_stream.save('n')


        elif exp == 'p-rabi':

            with program() as prog:

                a, n, I, Q = self.declare_vars([fixed, int, fixed, fixed])

                I_stream, Q_stream, n_stream = self.declare_streams(stream_num=3)

                update_frequency('qubit', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning) # sets the IF frequency of the qubit

                with for_(n, 0, n < n_avg, n + 1):
                    save(n, n_stream)
                    with for_(*from_array(a,var_arr)):  # Sweep pulse duration
                        play(pulse * amp(a), "qubit")
                        align("qubit", "rr")
                        measure("readout", "rr", None, *self.res_demod(I, Q))
                        save(I, I_stream)
                        save(Q, Q_stream)
                        wait(resettime_clk,"qubit")
                with stream_processing():
                    I_stream.buffer(len(var_arr)).average().save("I")
                    Q_stream.buffer(len(var_arr)).average().save("Q")
                    n_stream.save("n")

        elif exp == 'ramsey':

            with program() as prog:

               update_frequency('qubit', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning)

               n, t, I, Q = self.declare_vars([int, int, fixed, fixed])

               I_stream,Q_stream,n_stream = self.declare_streams(stream_num=3)

               with for_(n, 0, n < n_avg, n + 1):
                   save(n, n_stream)
                   with for_(*from_array(t,var_arr)):
                        play("pi_half", "qubit")
                        wait(t, "qubit")
                        # frame_rotation_2pi(phi, 'qubit')  # this was in Haimeng's code and was commented out by her,
                        play("pi_half", "qubit")
                        align("qubit","rr")
                        measure("readout", "rr", None, *self.res_demod(I, Q))
                        save(I, I_stream)
                        save(Q, Q_stream)
                        wait(resettime_clk, "qubit")

               with stream_processing():
                   I_stream.buffer(len(var_arr)).average().save("I")
                   Q_stream.buffer(len(var_arr)).average().save("Q")
                   n_stream.save("n")

        elif exp == 'echo':
            with program() as prog:

                n, t, I, Q = self.declare_vars([int, int, fixed, fixed])

                I_stream,Q_stream,n_stream = self.declare_streams(stream_num=3)

                update_frequency('qubit', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning)

                with for_(n, 0, n < n_avg, n + 1):
                    save(n,n_stream)
                    with for_(*from_array(t,var_arr)):
                        play("pi_half", "qubit")
                        wait(t, "qubit")
                        play("pi", "qubit")
                        wait(t, "qubit")
                        play("pi_half", "qubit")
                        align("qubit","rr")
                        measure("readout", "rr", None, *self.res_demod(I, Q))
                        save(I, I_stream)
                        save(Q, Q_stream)
                        wait(resettime_clk, "qubit")


                with stream_processing():
                    I_stream.buffer(len(var_arr)).average().save("I")
                    Q_stream.buffer(len(var_arr)).average().save("Q")
                    n_stream.save("n")

        elif exp == 'T1':
            with program() as prog:

                n, t, I, Q = self.declare_vars([int, int, fixed, fixed])

                I_stream, Q_stream, n_stream = self.declare_streams(stream_num=3)

                update_frequency('qubit', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning)

                with for_(n, 0, n < n_avg, n + 1):
                    save(n, n_stream)
                    with for_(*from_array(t,var_arr)):
                        play("pi", "qubit")
                        wait(t, 'rr')
                        align("qubit", "rr")
                        measure("readout", "rr", None,*self.res_demod(I,Q))
                        wait(resettime_clk, "qubit")
                        save(I, I_stream)
                        save(Q, Q_stream)

                with stream_processing():
                    I_stream.buffer(len(var_arr)).average().save("I")
                    Q_stream.buffer(len(var_arr)).average().save("Q")
                    n_stream.save('n')

        elif exp == 'ss':

            with program() as prog:

                i,n,I,Iexc,Q,Qexc = self.declare_vars([int,int,fixed,fixed,fixed,fixed])
                i_st, I_st,Q_st,I_st_exc,Q_st_exc = self.declare_streams(stream_num=5)

                with for_(i, 0, i < nIterations, i + 1):
                    with for_(n, 0, n < n_reps, n + 1):
                        # do nothing
                        align("qubit", "rr")
                        measure("readout", "rr", None, *self.res_demod(I, Q))
                        align('qubit','rr')
                        wait(resettime_clk, "qubit")
                        align('qubit','rr')
                        # apply pi-pulse
                        play("pi", "qubit")
                        # play("pi_half" * amp(0.01), "qubit", duration = clk(2e3))
                        # play("gauss"*amp(0.41/0.45), "qubit",duration=round(274*2/4))
                        align("qubit", "rr")
                        measure("readout", "rr", None, *self.res_demod(Iexc, Qexc))
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
                    I_st.buffer(n_reps).save('I')
                    Q_st.buffer(n_reps).save('Q')
                    I_st_exc.buffer(n_reps).save('Iexc')
                    Q_st_exc.buffer(n_reps).save('Qexc')
                    i_st.save('i')

        return prog

    #%%% pulse_exp
    def pulse_exp(self,sa = 0,
                      exp='rabi',
                      check_mixers=True,
                      n_avg = 2000,
                      tmin = 16,         # minimum pulse duration in nanoseconds
                      tmax = 10e3,    # maximum pulse duration in nanoseconds
                      dt = 500,        # step of sweep in nanoseconds
                      amp_q_scaling = 1,
                      fit = True,
                      plot = True,
                      detuning = 0e6):
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

        tmin = clk(tmin)
        tmax = clk(tmax)
        dt = clk(dt)
        t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)

        inst.set_attenuator(attenuation=self.pars['rr_atten'])
        # self.update_value('qubit_LO', self.pars['qubit_LO'])
        # self.update_value('rr_LO', self.pars['rr_LO'])

        # qb_IF = self.pars['qubit_freq']-self.pars['qubit_LO']
        # self.update_value('qubit_IF', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning)
        self.check_mix_cal(sa, amp_q = amp_q_scaling,check = check_mixers, threshold = -55)

        prog = self.make_sequence(exp=exp,var_arr=t_arr,detuning=detuning,n_avg=n_avg,amp_q_scaling=amp_q_scaling)

        datadict, job = self.get_results(jobtype = prog, result_names = ["I", "Q"], n_total=n_avg, notify = False)

        qb_power = self.get_power(sa,freq=self.pars['qubit_LO']+self.pars['qubit_IF'],reference=0,amp_q = amp_q_scaling, span=1e6,config=True,output=False)

        # t_arr = np.array(datadict["t"])/1e3 # times in microseconds
        if exp == 'echo':
            t_arr = np.array(t_arr)*4/1e3 * 2
        else:
            t_arr = np.array(t_arr)*4/1e3
        I = np.array(datadict["I"])
        Q = np.array(datadict["Q"])
        ydata = np.abs(I+1j*Q)

        if fit:
            fitted_pars, error = pf.fit_data(t_arr,ydata,sequence=exp,dt=t_arr[-1]*1e-6/len(t_arr))
            if plot:
                pf.plot_data(t_arr,ydata,sequence=exp,fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=self.pars['pi_half_len'],
                             qubitDriveFreq=self.pars['qubit_LO']+self.pars['qubit_IF'],qb_power = qb_power,iteration=iteration)
        elif plot:
            fitted_pars = None
            pf.plot_data(t_arr, ydata, sequence = exp)
        else:
            fitted_pars = None




        if exp == 'rabi':
            self.update_value('pi_half_len',4*(clk(fitted_pars[1]/4)+1))

        # self.update_value('qubit_IF',qb_IF)

        exp_dict = {'date/time':     datetime.now(),
                   'nAverages': n_avg,
                         'Tmax': tmax,
                         'dt':   dt,
                         'pi2': self.pars['pi_half_len'],
                         'A_d':     amp_q_scaling,
                         'w_d':  self.pars['qubit_LO']+self.pars['qubit_IF']-detuning,
                         'w_LO': self.pars['qubit_LO'],
                'wait_period':  self.pars['qubit_resettime'],
                }

        # save data
        with open(f"D:\weak_measurements\{exp}\data_{iteration:03d}.csv","w") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(exp_dict.keys())
            writer.writerow(exp_dict.values())
            writer.writerow(t_arr)
            writer.writerow(I)
            writer.writerow(Q)

        return t_arr, I, Q, job, fitted_pars





#%% CALIBRATIONS
  #%%% tof_cal
    def tof_cal(self):
        qmm = QuantumMachinesManager(host=host, port=port)

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

#%%% cal_pi_pulse
    def cal_pulse(self, pulse = 'pi', pulse_len_target = 20, starting_amp = 0.44, **rabi_kwargs):

        if pulse == 'pi':
            scaling = 2
        elif pulse == 'pi_half':
            scaling = 4

        print(f"CALIBRATING {pulse} pulse... target time = {pulse_len_target}")

        # start with a clean slate based on input parameters
        self.update_value(f'{pulse}_amp', starting_amp)
        self.update_value(f'{pulse}_len', pulse_len_target)

        # run power rabi and extract fit parameter A, which is the scaling factor
        print("Calibrating fit...")
        amps, I, Q, job, A = self.power_rabi(check_mixers = False, pulse = pulse, fit = True, plot = True, **rabi_kwargs)

        print(f"Fitted scaling factor A/{scaling} = {round(A/scaling, ndigits=4)}")

        # check: if A / 2 < 1, then you've probably found your pi pulse and you're donee
        # if A / 2 >= 1, then you probably need to increrase your pi_half_len_target
        # this code will automatically recalibrate with a new pi_half_len_target based on the estimate of A
        if A / scaling < 1:

            print(f"{pulse} pulse calibration successful, updating config...")

            # update pi_amp
            self.update_value(f'{pulse}_amp', starting_amp * A / scaling)

            # as a check to your pi pulse, run single shot and run power rabi again to see your pi pulse in action
            if pulse == 'pi':
                self.single_shot(nIterations=1)
            amps, I, Q, job, A = self.power_rabi(check_mixers = False, pulse = pulse, fit = True, plot = True, **rabi_kwargs)

        else:

            # your new pi pulse time to try is the old time, scaled by A/2, times 1.1 to overshoot a bit
            # convert to clock cycles then multiply by 4 to ensure a pulse time w/ integer clock cycles
            newtarget = 4 * clk(pulse_len_target * (A/scaling) * 1.1)

            print(f"The {pulse} length target time was too short. Trying again for target time {newtarget} ns.")

            self.cal_pulse(pulse_len_target = newtarget, pulse = pulse, starting_amp = starting_amp, **rabi_kwargs)


    #%%% cal_freq
    def cal_freq(self,  min_delta=10e3,
                        check_mixers = False,
                        recalibrate_pulses = True,
                        dt_scaling = 4,
                        **kwargs):

        print("RUNNING initial Ramsey experiment...")
        t_arr, I, Q, job, fitted_pars = self.pulse_exp(exp='ramsey', check_mixers = check_mixers,dt=40,tmin=20,tmax=2e4,detuning=0,**kwargs)

        delta = abs(fitted_pars[1]) * 1e6
        tau = fitted_pars[3] * 1e3
        tmax = 2 * tau

        self.update_value('qubit_freq', self.pars['qubit_freq'] + delta)

        n = 0

        while n < 10:

            n += 1

            print(f'RUNNING Ramsey experiment, delta = {round(delta/1e3, ndigits=6)} kHz')

            if delta < 100e3:
                scaling = dt_scaling * 8
            else:
                scaling = dt_scaling

            dt = round(get_dt(target_fs=delta)/scaling)

            t_arr, I, Q, job, fitted_pars = self.pulse_exp(exp='ramsey', detuning = 0e6, check_mixers = check_mixers,dt=dt,tmin=20,tmax=tmax, **kwargs)
            delta_result = abs(fitted_pars[1]) * 1e6

            print(f'RESULT: delta_result = {round(delta_result/1e3, ndigits=6)} kHz')

            if delta_result < min_delta:
                print(f'SUCCESS: delta_result < min_delta = {round(min_delta/1e3, ndigits=6)} kHz :) ')
                rounded_freq = round_from_uncertainty(self.pars['qubit_freq'], 1e3, scale=1e9) # round frequency to nearest kHz
                self.update_value('qubit_freq', rounded_freq)

                if recalibrate_pulses:
                    self.cal_pulse(pulse = 'pi', pulse_len_target = self.pars['pi_len'], **kwargs)
                    self.cal_pulse(pulse = 'pi_half', pulse_len_target = self.pars['pi_half_len'], **kwargs)

                break
            elif delta_result > delta:
                print(f'WRONG DIRECTION: delta_result > delta')
                self.update_value('qubit_freq', self.pars['qubit_freq'] - 2 * delta)
                delta = delta_result
            else:
                print(f'RIGHT DIRECTION: delta_result < delta')
                self.update_value('qubit_freq', self.pars['qubit_freq'] + delta_result)
                delta = delta_result



        print(f"Didn't converge within min_freq = {minfreq}")



    #%%% cal_freq_new
    def cal_freq_new(self,  min_delta=10e3,
                        check_mixers = False,
                        recalibrate_pulses = True,
                        dt_scaling = 4,
                        num_tries = 10,
                        **kwargs):

        # initial Ramsey experiment to find current detuning
        print("RUNNING initial Ramsey experiment...")
        t_arr, I, Q, job, fitted_pars = self.pulse_exp(exp='ramsey', check_mixers = check_mixers,dt=40,tmin=20,tmax=2e4,detuning=0,**kwargs)

        delta = abs(fitted_pars[1]) * 1e6
        tau = fitted_pars[3] * 1e3
        tmax = 2 * tau

        delta_dict = { delta : {"freq" : self.pars['qubit_freq'], "tried_increase" : False}}
        #


        for n in np.arange(0, num_tries - 1):

            # update value of qubit_freq to be that corresponding to the smallest delta in delta_dict,
            # plus or minus delta
            delta = min(delta_dict.keys())
            old_freq = delta_dict[delta]["freq"]
            tried_increase = delta_dict[delta]["tried_increase"]
            if not tried_increase:
                self.update_value('qubit_freq', old_freq + delta)
                delta_dict[delta]["tried_increase"] = True
            else:
                self.update_value('qubit_freq', old_freq - delta)

            # choose appropriate dt for the frequency we're trying to resolve
            if delta < 100e3:
                scaling = dt_scaling * 8
            else:
                scaling = dt_scaling

            dt = round(get_dt(target_fs=delta)/scaling)

            # run Ramsey experiment and get new delta
            t_arr, I, Q, job, fitted_pars = self.pulse_exp(exp='ramsey', detuning = 0e6, check_mixers = check_mixers,dt=dt,tmin=20,tmax=tmax, **kwargs)
            delta_result = abs(fitted_pars[1]) * 1e6

            if delta_result < min_delta:
                print(f'SUCCESS: delta_result < min_delta = {round(min_delta/1e3, ndigits=6)} kHz :) ')
                rounded_freq = round_from_uncertainty(self.pars['qubit_freq'], 1e2, scale=1e9) # round frequency to nearest kHz
                self.update_value('qubit_freq', rounded_freq)

                if recalibrate_pulses:
                    self.cal_pulse(pulse = 'pi', pulse_len_target = self.pars['pi_len'], **kwargs)
                    self.cal_pulse(pulse = 'pi_half', pulse_len_target = self.pars['pi_half_len'], **kwargs)

                return # if we were successful, we skip over the statement at the end
            # ONLY if we did better than before, add the new delta_result and qubit freq to the dict
            elif delta_result < delta:
                delta_dict[delta_result] = {"freq" : self.pars["qubit_freq"], "tried_increase" : False}

        delta = min(delta_dict.keys())
        freq = delta_dict[delta]["freq"]
        self.update_value('qubit_freq', freq)

        print(f"Didn't converge within min_freq = {minfreq} within {num_tries} tries. The best result was {round(delta/1e3, digits=3)} kHz at qubit freq {round(freq / 1e9, 7)} GHz.")


#%%% calibrate
    def calibrate(self,   pulse_len_target = 48,
                          starting_amp = 0.44,
                          min_delta = 10e3,
                          recalibrate_pulses = True,
                          T1_kwargs     = {"dt": 1e3, "tmin" : 100, "tmax" : 100e3},
                          ramsey_kwargs = {"dt": 200, "tmin" : 20, "tmax" : 40e3},
                          echo_kwargs   = {"dt": 1e3, "tmin" : 50, "tmax" : 120e3},
                          **kwargs
                          ):

        # time of flight
        self.tof_cal() # this can be modified to actually calibrate tof automatically, rather than just display plot

        # pi pulse length / time
        self.cal_pulse(pulse = 'pi', pulse_len_target = pulse_len_target, starting_amp = starting_amp, **kwargs)
        self.cal_pulse(pulse = 'pi_half', pulse_len_target = pulse_len_target, starting_amp = starting_amp, **kwargs)

        # frequency
        self.cal_freq(min_delta = min_delta, recalibrate_pulses = recalibrate_pulses, **kwargs)

        # T1, T2 ramsey, T2 echo
        print("RUNNING longer fit experiments...")
        t_arr, I, Q, job, fitted_pars = self.pulse_exp(exp='T1', check_mixers=True, **T1_kwargs, **kwargs)
        T1 = fitted_pars[1]

        t_arr, I, Q, job, fitted_pars = self.pulse_exp(exp='ramsey', check_mixers=True, **ramsey_kwargs, **kwargs)
        T2_ramsey = fitted_pars[3]

        t_arr, I, Q, job, fitted_pars = self.pulse_exp(exp='echo', check_mixers=True, **echo_kwargs, **kwargs)
        T2_echo = fitted_pars[1]

        print(f"CALIBRATION SUCCESSFUL. Qubit parameters: \nT1 = {T1}, \nT2 ramsey = {T2_ramsey}, \nT2_echo = {T2_echo} ")





#%% HELPERS
#%%% Get Results
#%%%% get_results
    def get_results(self,jobtype, result_names = ["I", "Q"],
                                    progress_key = "n",
                                    showprogress = True,
                                    n_total = 1000,
                                    notify = False,
                                    liveplot = False):
        """

        Args:
            jobtype (QUA program): QUA program to execute.
            result_names (array of strings, optional): The names of the variables we want to stream. Defaults to ["I", "Q", "n"].
            showprogress (TYPE, optional): DESCRIPTION. Defaults to False.
            n_total (TYPE, optional): DESCRIPTION. Defaults to 1000.
            notify (TYPE, optional): DESCRIPTION. Defaults to False.
            liveplot (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            datadict (TYPE): DESCRIPTION.
            TYPE: DESCRIPTION.

        """


        # Open Communication with the Server
        qmm = QuantumMachinesManager(host=host, port=port)

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

            I_handle = handles_dict["I"]

            for handle in handles_dict.values():
                handle.wait_for_values(1)
                is_processing = lambda: handle.is_processing()
            i0 = 0
            while(I_handle.is_processing()): # while(is_processing()):
                i = handles_dict["i"].fetch_all() # retrieve iteration value n
                datadict = self.get_data_from_handles(handles_dict)
                pf.plot_single_shot(datadict,axes=ax)

        if showprogress:
            n_handle = res_handles.get(progress_key)
            make_progress_meter(n_handle, n_total)
#
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

    #%%%% get_data_from_handles
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

#%%% Mixer
    #%%%% config_sa
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

    #%%%% get_power
    def get_power(self,sa,freq,reference=-100,span=1e6, amp_q=1, config=False, plot=False, output = True):
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
        inst.set_attenuator(attenuation=0)

        # skips configuring the spectrum analyzer. Used only when optimizing mixer
        if config:
            self.play_pulses(amp_q=amp_q)
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

        inst.set_attenuator(attenuation=self.pars['rr_atten'])

        return power

    #%%%% opt_mixer
    def opt_mixer(self,sa,cal,mode,element,freq_span=1e6,amp_q=1,reference = -30, plot=True):
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
        qm = self.play_pulses(amp_q=amp_q)

        if cal == 'LO':
            freq = self.pars[f'{element}_LO']
            par1 = self.pars[f'{element}_mixer_offsets'][0]
            par2 = self.pars[f'{element}_mixer_offsets'][1]
            print(f'LO at {round(freq*1e-9,5)} GHz\nCurrent I_offset = {round(par1*1e3,1)} mV, Current Q_offset = {round(par2*1e3,1)} mV')
        elif cal == 'SB':
            freq = self.pars[f'{element}_LO'] - self.pars[f'{element}_IF']
            par1 = self.pars[f'{element}_mixer_imbalance'][0]
            par2 = self.pars[f'{element}_mixer_imbalance'][1]
            # if np.abs(par2) > 150e-3:
            #     par2 = 0
            print(f'Sideband at {round(freq*1e-9,5)} GHz\nCurrent gain = {round(par1,4)}, Current phase = {round(par2,4)}')

        if element == 'rr':
            inst.set_attenuator(attenuation=0)

        self.config_sa(sa,freq,span=freq_span,reference=reference) # setup spectrum analyzer for measurement

        # initialize sweep parameters
        if cal == 'LO':
            if mode == 'coarse':
                span = 40e-3
                step = 10e-3
            elif mode == 'intermediate':
                span = 5e-3
                step = 1e-3
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
                span = 6e-3
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
                        qm.set_mixer_correction(element,int(self.pars[f'{element}_IF']), int(self.pars[f'{element}_LO']), self.IQ_imbalance(par1, par2))
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
            qm.set_mixer_correction(element,int(self.pars[f'{element}_IF']), int(self.pars[f'{element}_LO']), self.IQ_imbalance(par1_arr[argmin[0]],par2_arr[argmin[1]]))
            self.update_value(f'{element}_mixer_imbalance',[par1_arr[argmin[0]],par2_arr[argmin[1]]])
            print(f"optimal gain = {round(par1_arr[argmin[0]],4)}, optimal phi = {round(par2_arr[argmin[1]],4)}")

        if element == 'rr':
            inst.set_attenuator(self.pars['rr_atten'])

        print(f'Power: {np.amin(power_data)} dBm at {freq/1e9} GHz')
        if plot:
            pf.plot_mixer_opt(par1_arr, par2_arr, power_data,cal=cal,element=element,fc=freq)

    #%%%% check_mixer_calibration
    def check_mix_cal(self, sa, amp_q = 1, check = True, threshold = -50):
        # checks LO leakage and image sideband power and opt mixer accordingly
        if check:
            for element in ['qubit','rr']:
                for j in ['LO','SB']:
                    if j == 'LO':
                        fc = self.pars[f'{element}_LO']
                    elif j == 'SB':
                        fc = self.pars[f'{element}_LO'] -  self.pars[f'{element}_IF']
                    print(f'Checking {element} {j}')
                    leak = self.get_power(sa, freq = fc, reference = -20, config = True, amp_q = amp_q, plot = False)
                    while leak > threshold:
                        leak0 = leak
                        ref = leak + 10
                        print(f'Minimizing {element} {j} leakage')
                        if leak > - 40:
                            self.opt_mixer(sa,cal = j, mode = 'coarse', amp_q = amp_q, element = element, reference = ref)
                        elif leak < - 40 and leak > - 55:
                            self.opt_mixer(sa, cal = j, mode = 'intermediate', amp_q = amp_q, element = element, reference = ref)
                        elif leak < - 55:
                            self.opt_mixer(sa, cal = j, mode = 'fine', amp_q = amp_q, element = element, reference = ref)

                        leak = self.get_power(sa,freq = fc,reference = ref, amp_q = amp_q, config = True, plot = False)

                        if np.abs(leak - leak0) < 5:
                            print("Can't optimize mixer further")
                            break
        else:
            pass

#%%% Macros

    #%%%% declare_vars
    def declare_vars(self,types):
        return [declare(tp) for tp in types]

    #%%%% declare_streams
    def declare_streams(self,stream_num=1):
        return [declare_stream() for num in range(stream_num)]

    #%%%% res_demod
    def res_demod(self,I, Q):

        return (dual_demod.full("integW_cos", "out1", "integW_minus_sin", "out2", I),
                dual_demod.full("integW_sin", "out1", "integW_cos", "out2", Q))







    #%% CONFIGURATION

    #%%% Helpers

    #%%%% IQ_imbalance
    def IQ_imbalance(self,g, phi):
        c = np.cos(phi)
        s = np.sin(phi)
        N = 1 / ((1-g**2)*(2*c**2-1))
        return [float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]

    #%%%% delayed_gauss
    def delayed_gauss(self,amp=0.45, length=4, sigma=2):
        gauss_arg = np.linspace(-sigma, sigma, length)
        delay = 16 - length - 4
        if delay < 0:
            return amp * np.exp(-(gauss_arg ** 2) / 2)

        return np.r_[np.zeros(delay), amp * np.exp(-(gauss_arg ** 2) / 2), np.zeros(4)]

    #%%%%
    def init_instruments(self):

        inst.set_qb_LO(self.pars['qubit_LO'])
        inst.set_rr_LO(self.pars['rr_LO'])
        inst.set_attenuator(self.pars['rr_atten'])


    #%%%% update_value
    def update_value(self,key,value):
        print(f'Updating {key} to {value}')
        self.pars[key] = value
        self.make_config(self.pars)

        if key == 'qubit_LO':
            inst.set_qb_LO(value)
        elif key == 'rr_LO':
            inst.set_rr_LO(value)
        elif key == 'rr_atten':
            inst.set_attenuator(value)
        elif key == 'qubit_freq':
            self.update_value('qubit_IF',self.pars['qubit_freq']-self.pars['qubit_LO'])

        self.write_pars()

    #%%%% write_pars
    def write_pars(self):
        with open(f'{self.name}_pars.json', "w") as outfile:
            json.dump(self.pars, outfile)

    #%%%% remove_key
    def remove_key(self, key):
        print(f'Removing {key} from pars')
        del self.pars[key]
        self.write_pars()

    #%%%% add_key
    def add_key(self, key, value):
        print(f'Adding {key} = {value} to pars')
        self.pars[key] = value
        self.write_pars()

    #%%% Make Config
    #%%%% make_config
    def make_config(self, pars):
        gauss_wf_4ns = self.delayed_gauss()

        self.config = {

            "version": 1,

            "controllers": {
                "con2": {
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
                    "intermediate_frequency": pars['qubit_IF'],
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
                    "intermediate_frequency": pars['rr_IF'],
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
                    "length": pars['pi_len'],
                    "waveforms": {
                        "I": "pi_wf_i1",
                        "Q": "pi_wf_q1",
                    },
                },
                "Xpi_pulse": {
                    "operation": "control",
                    "length": pars['pi_len'],
                    "waveforms": {
                        "I": "pi_wf_i1",
                        "Q": "pi_wf_q1",
                    },
                },
                "Ypi_pulse": {
                    "operation": "control",
                    "length": pars['pi_len'],
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
                "pi_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_amp'] * gaussian(pars['pi_len'], pars['pi_len']/5)]},
                "pi_wf_q1": {"type": "constant", "sample": 0.0},
                "pi_half_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_half_amp'] * gaussian(pars['pi_half_len'], pars['pi_half_len']/5)]},
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
                "qubit": [{"intermediate_frequency": pars['qubit_IF'], "lo_frequency": pars['qubit_LO'], "correction": self.IQ_imbalance(*pars['qubit_mixer_imbalance'])}],
                "rr": [{"intermediate_frequency": pars['rr_IF'], "lo_frequency": pars['rr_LO'], "correction": self.IQ_imbalance(*pars['rr_mixer_imbalance'])}],
            }
        }

#%% GRAVEYARD

#%%% MACROS
    #%%%% pi_pulse
    # play pi pulse, measure, and wait
    def pi_pulse(self, I, Q, resettime = 400e3):

        align("qubit", "rr")
        play("pi", "qubit")
        align("qubit", "rr")
        measure("readout", "rr", None, *self.res_demod(I, Q))
        wait(clk(resettime), "qubit")
        align("qubit", "rr")

        return I, Q

    #%%%% save_many
    # this needs to be tested...
    # save multiple variables to multiple streams in one line

    # example: instead of writing
    # save(I, I_st)
    # save(Q, Q_st)

    # can write
    # self.save_many([I, I_st], [Q, Q_st])
    def save_many(*pairs):

        for (var, var_stream) in pairs:

            save(var, var_stream)


 #%%% UNIT TESTS
    #%%%% pi pulse
    def test_pi_pulse(self, n_avg = 5000, tolerance = 1e-3):

        with program() as PiTest:

            n, I, Q = self.declare_vars([int, fixed, fixed])
            I_st, Q_st = self.declare_streams(stream_num=2)

            with for_(n, 0, n < n_avg, n + 1):
                wait(clk(100e3),'rr')
                # first, measure without doing pi pulse
                measure("readout", "rr", None, *self.res_demod(I, Q))
                wait(clk(100e3), "rr")
                save(I, I_st)
                save(Q, Q_st)

                # then, measure after doing pi pulse
                I, Q = self.pi_pulse(I, Q)
                save(I, I_st)
                save(Q, Q_st)

            with stream_processing():
                I_st.buffer(2).average().save("I")
                Q_st.buffer(2).average().save("Q")

        # run job
        datadict, _ = self.get_results(jobtype = PiTest,
                                         result_names = ["I", "Q"],
                                         n_total=n_avg,
                                         showprogress=False,
                                         notify = False)
        # process data from job
        I0, I1 = datadict["I"]
        Q0, Q1 = datadict["Q"]
        p0 = I0 + 1j * Q0
        p1 = I1 + 1j * Q1
        mag1 = round(np.abs(p1) * 1e3, ndigits=3)
        mag0 = round(np.abs(p0) * 1e3, ndigits=3)
        dist = round(mag1 - mag0, ndigits=3)

        print(f"Distinguishability: {mag1} mV - {mag0} mV = {dist} mV")

        ### eventually, include code that checks that this matches what
        # we expect from power rabi to within {tolerance}

        return datadict

#%%%% test helpers

    # e.g. operation = "pi"
    def get_pulse(self, operation, element):
        return self.config["elements"][element]["operations"][operation]

    def get_waveform(self, operation, element, quad = "I"):
        pulse = self.get_pulse(operation, element)
        return self.config["pulses"][pulse]["waveforms"][quad]

    def waveform(self, operation, element):
        waveform_name = self.get_waveform(operation, element)
        return self.config["waveforms"][waveform_name]["samples"]


    #%%% IQ_blobs
    """
IQ_blobs.py: Measure the qubit in the ground and excited state to create the IQ blobs.
If the separation and the fidelity are good enough, gives the parameters needed for active reset
"""


# ###################
# # The QUA program #
# ###################

# n_runs = 10000

# cooldown_time = 5 * qubit_T1 // 4

# with program() as IQ_blobs:
#     n = declare(int)
#     I_g = declare(fixed)
#     Q_g = declare(fixed)
#     I_g_st = declare_stream()
#     Q_g_st = declare_stream()
#     I_e = declare(fixed)
#     Q_e = declare(fixed)
#     I_e_st = declare_stream()
#     Q_e_st = declare_stream()

#     with for_(n, 0, n < n_runs, n + 1):
#         measure(
#             "readout",
#             "resonator",
#             None,
#             dual_demod.full("rotated_cos", "out1", "rotated_sin", "out2", I_g),
#             dual_demod.full("rotated_minus_sin", "out1", "rotated_cos", "out2", Q_g),
#         )
#         save(I_g, I_g_st)
#         save(Q_g, Q_g_st)
#         wait(cooldown_time, "resonator")

#         align()  # global align

#         play("pi", "qubit")
#         align("qubit", "resonator")
#         measure(
#             "readout",
#             "resonator",
#             None,
#             dual_demod.full("rotated_cos", "out1", "rotated_sin", "out2", I_e),
#             dual_demod.full("rotated_minus_sin", "out1", "rotated_cos", "out2", Q_e),
#         )
#         save(I_e, I_e_st)
#         save(Q_e, Q_e_st)
#         wait(cooldown_time, "resonator")

#     with stream_processing():
#         I_g_st.save_all("I_g")
#         Q_g_st.save_all("Q_g")
#         I_e_st.save_all("I_e")
#         Q_e_st.save_all("Q_e")

# #####################################
# #  Open Communication with the QOP  #
# #####################################
# qmm = QuantumMachinesManager(qop_ip)

# qm = qmm.open_qm(config)

# job = qm.execute(IQ_blobs)
# res_handles = job.result_handles
# res_handles.wait_for_all_values()
# Ig = res_handles.get("I_g").fetch_all()["value"]
# Qg = res_handles.get("Q_g").fetch_all()["value"]
# Ie = res_handles.get("I_e").fetch_all()["value"]
# Qe = res_handles.get("Q_e").fetch_all()["value"]

# angle, threshold, fidelity, gg, ge, eg, ee = two_state_discriminator(Ig, Qg, Ie, Qe, b_print=True, b_plot=True)

#########################################
# The two_state_discriminator gives us the rotation angle which makes it such that all of the information will be in
# the I axis. This is being done by setting the `rotation_angle` parameter in the configuration.
# See this for more information: https://qm-docs.qualang.io/guides/demod#rotating-the-iq-plane
# Once we do this, we can perform active reset using:
#########################################
#
# # Active reset:
# with if_(I > threshold):
#     play("pi", "qubit")
#
#########################################
#
# # Active reset (faster):
# play("pi", "qubit", condition=I > threshold)
#
#########################################
#
# # Repeat until success active reset
# with while_(I > threshold):
#     play("pi", "qubit")
#     align("qubit", "resonator")
#     measure("readout", "resonator", None,
#                 dual_demod.full("rotated_cos", "out1", "rotated_sin", "out2", I))
#
#########################################
#
# # Repeat until success active reset, up to 3 iterations
# count = declare(int)
# assign(count, 0)
# cont_condition = declare(bool)
# assign(cont_condition, ((I > threshold) & (count < 3)))
# with while_(cont_condition):
#     play("pi", "qubit")
#     align("qubit", "resonator")
#     measure("readout", "resonator", None,
#                 dual_demod.full("rotated_cos", "out1", "rotated_sin", "out2", I))
#     assign(count, count + 1)
#     assign(cont_condition, ((I > threshold) & (count < 3)))
#
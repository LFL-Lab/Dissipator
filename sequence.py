# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:00:18 2023

@author: lfl
"""
from qm.qua import *
from qualang_tools.loops import from_array
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm import SimulationConfig
import numpy as np
from Utilities import clk
from helper_functions import res_demod, declare_vars, declare_streams


class sequence():
    
    # def __init__(self, name, n_avg=100, amp_r_scale=1, amp_ffl_scale=1, **kwargs):
    def __init__(self,qb,name,**kwargs):
        self.name = name
        # self.seq_pars = qb.seq_pars
        self.qb_pars = qb.pars
        self.seq_pars = {}
        # self.exp_dict = {'n_avg': n_avg,
        #                  'amp_r_scale': amp_r_scale,
        #                  'amp_ffl_scale': amp_ffl_scale}
        for key in kwargs.keys():
            self.seq_pars[key] = kwargs.get(key)
        
        
    # def make_sequence(self, qb, tmin = 16, tmax=1e3, dt=8, scrambling_amp=0., ffl_len=2e3, saturation_dur=2e3, with_ffl='True', with_scram='True', detuning=0):
        

    def single_tone_spectroscopy(self):
        """Single tone spectroscopy sequence for linear resonators"""
        n_avg = self.qb_pars['n_avg']
        IF_min = self.seq_pars['IF_min']
        IF_max = self.seq_pars['IF_max']
        df = self.seq_pars['df']
        freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
        res_ringdown_time = self.qb_pars['resettime']['rr']
        tof_in_us = self.qb_pars['tof'] * 1e-3
        IQ_rotation = self.qb_pars['IQ_rotation']
        phase_prefactor = tof_in_us * 0.000064 # 0.000064 = 2^6 * 1e-6, to compensate for right-bitshift by 6 bits
        initial_frame_rotation = self.qb_pars['tof'] * IF_min/1e9
        ### QUA code ###
        with program() as prog:
            n, I, Q, f = declare_vars([int, fixed, fixed, int])
            I_st, Q_st, n_stream = declare_streams(stream_num=3)
            phase_shift = declare(fixed)
            df = declare(int, value=int(df))

            with for_(n, 0, n < n_avg, n + 1):
                save(n,n_stream)
                reset_frame('rr')
                frame_rotation_2pi(IQ_rotation + initial_frame_rotation,'rr')
                # with for_each_(f, freqs_list):
                with for_(*from_array(f, freqs)):
                    update_frequency("rr", f,keep_phase=False)
                    assign(phase_shift, Cast.mul_fixed_by_int(phase_prefactor, (df>>6)))
                    frame_rotation_2pi(phase_shift,'rr')
                    wait(res_ringdown_time, "rr")
                    measure("readout", "rr", None,*res_demod(I, Q,switch_weights=self.qb_pars['switch_weights']))
                    save(I, I_st)
                    save(Q, Q_st)


            with stream_processing():
                I_st.buffer(len(freqs)).average().save('I')
                Q_st.buffer(len(freqs)).average().save('Q')
                n_stream.save('n')

        return prog
                    
    def two_tone_spectroscopy(self,):
        """Two tone resonator spectroscopy with power sweep. Can be used for qubit or photon-resolved resonator spectroscopy"""
        n_avg = self.qb_pars['n_avg']
        freqs = np.arange(self.seq_pars['IF_min'], self.seq_pars['IF_max'] + self.seq_pars['df']/2, self.seq_pars['df'], dtype=int)
        target_elem = self.seq_pars['target_res'] # the resonator to be probed, choices are 'qubit', 'rr' or 'cavity'
        readout_elem = self.seq_pars['readout_res'] # the resonator to be readout, choices are 'rr' or 'cavity' or 'qubit'
        
        print(f"Performing two-tone pulsed spectroscopy on {target_elem} using {readout_elem} for readout")
        a = self.seq_pars['amp_q_scaling']
        saturation_dur = self.seq_pars['saturation_dur']
        
        on_off = self.seq_pars['on_off']
        # amp_ffl_scale = self.exp_dict['amp_ffl_scale']
        reset_time_target = self.qb_pars['resettime'][target_elem]
        reset_time_readout = self.qb_pars['resettime'][readout_elem]
        print(f"Reset time for {target_elem} is {reset_time_target} and for {readout_elem} is {reset_time_readout}")

        with program() as prog:
            # update_frequency('rr', self.qb_pars['rr_IF'])
            n, I, Q, f = declare_vars([int, fixed, fixed, int])
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            if on_off:
                I_b, Q_b, I_tot,Q_tot = declare_vars([fixed,fixed, fixed, fixed])
            
        
            with for_(n, 0, n < n_avg, n + 1):
                save(n, n_stream)
                with for_(*from_array(f, freqs)):
                    update_frequency(target_elem, f)
                    # measure background
                    if on_off:
                        measure("readout", readout_elem, None, *res_demod(I_b, Q_b,switch_weights=self.qb_pars['switch_weights']))
                        align(readout_elem, target_elem) # wait for operations on readout resonator to finish before playing qubit pulse
                        wait(reset_time_readout, readout_elem)
                        align(readout_elem, target_elem)

                    play('const'*amp(a), target_elem, duration=saturation_dur)
                    align(target_elem, readout_elem)
                    measure("readout", readout_elem, None, *res_demod(I, Q,switch_weights=self.qb_pars['switch_weights']))
                    wait(reset_time_target,target_elem)
             
                    if on_off:
                        assign(I_tot, I - I_b) 
                        assign(Q_tot, Q - Q_b)
                        save(I_tot, I_stream)
                        save(Q_tot, Q_stream)
                    else:
                        save(I, I_stream)
                        save(Q, Q_stream)
        
            with stream_processing():
                I_stream.buffer(len(freqs)).average().save('I')
                Q_stream.buffer(len(freqs)).average().save('Q')
                n_stream.save('n')

        return prog
    
    def two_tone_spectroscopy_amp_sweep(self,):
        """Two tone resonator spectroscopy with power sweep. Can be used for qubit or photon-resolved resonator spectroscopy"""
        n_avg = self.qb_pars['n_avg']
        freqs = np.arange(self.seq_pars['IF_min'], self.seq_pars['IF_max'] + self.seq_pars['df']/2, self.seq_pars['df'], dtype=int)
        target_elem = self.seq_pars['target_res'] # the resonator to be probed, choices are 'qubit', 'rr' or 'cavity'
        readout_elem = self.seq_pars['readout_res'] # the resonator to be readout, choices are 'rr' or 'cavity' or 'qubit'
        amp_arr = np.arange(self.seq_pars['amin'], self.seq_pars['amax'] + self.seq_pars['da']/2, self.seq_pars['da'], dtype=float)
        print(amp_arr)
        qubit_amp = self.seq_pars['amp_q_scaling']
        saturation_dur = self.seq_pars['saturation_dur']
        on_off = self.seq_pars['on_off']
        # amp_ffl_scale = self.exp_dict['amp_ffl_scale']
        reset_time = self.qb_pars['resettime'][target_elem]

        with program() as prog:
            # update_frequency('rr', self.qb_pars['rr_IF'])
            n, I, Q = declare_vars([int, fixed, fixed])
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            if on_off:
                I_b, Q_b, I_tot,Q_tot = declare_vars([fixed,fixed, fixed, fixed])
            f,a = declare_vars([int,fixed])

            with for_(*from_array(a,amp_arr)):
                with for_(n, 0, n < n_avg, n + 1):
                    save(n, n_stream)
                    with for_(*from_array(f, freqs)):
                        update_frequency(target_elem, f)
                        # measure background
                        if on_off:
                            measure("readout", readout_elem, None, *res_demod(I_b, Q_b,switch_weights=self.qb_pars['switch_weights']))
                            wait(reset_time, target_elem)
                            align(readout_elem, target_elem) # wait for operations on readout resonator to finish before playing qubit pulse

                        play('const'*amp(a), target_elem, duration=saturation_dur)
                        align(target_elem, readout_elem)
                        measure("readout", readout_elem, None, *res_demod(I, Q,switch_weights=self.qb_pars['switch_weights']))
                        wait(reset_time,target_elem)

                        if on_off:
                            assign(I_tot, I - I_b) 
                            assign(Q_tot, Q - Q_b)
                            save(I_tot, I_stream)
                            save(Q_tot, Q_stream)
                        else:
                            save(I, I_stream)
                            save(Q, Q_stream)
            
                with stream_processing():
                    I_stream.buffer(len(freqs)).average().save('I')
                    Q_stream.buffer(len(freqs)).average().save('Q')
                    n_stream.save('n')

        return prog
    
    def power_rabi(self,):
        amin = self.seq_pars['amin']
        amax = self.seq_pars['amax']
        da = self.seq_pars['da']
        amp_arr = np.arange(amin, amax + da/2, da, dtype=float)
        resettime_clk= self.qb_pars['resettime']['qubit']
        n_avg = self.qb_pars['n_avg']
        amp_r_scale = self.qb_pars['amp_r_scale']
        # amp_ffl_scale = self.qb_pars['amp_ffl_scale']
        rr_IF = self.qb_pars['rr_IF']
        # ffl_IF = self.qb_pars['ffl_IF']
        with program() as prog:
            update_frequency('rr', rr_IF) 
            # update_frequency('ffl', ffl_IF) 
            n, a, I, Q = declare_vars([int, fixed, fixed, fixed])
        
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            
            with for_(n, 0, n < n_avg, n + 1):
                save(n, n_stream)
                with for_(*from_array(a, amp_arr)):
                    play("gauss"*amp(a), "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights'])) 
                    wait(resettime_clk, "qubit")
                    save(I, I_stream)
                    save(Q, Q_stream)
        
            with stream_processing():
                I_stream.buffer(len(amp_arr)).average().save("I")
                Q_stream.buffer(len(amp_arr)).average().save("Q")
                n_stream.save('n')   

        return prog
    

    def ramsey(self,):
        tmin = self.seq_pars['tmin']
        tmax = self.seq_pars['tmax']
        dt = self.seq_pars['dt']
        t_arr = np.arange(tmin, tmax + dt/2, dt, dtype=float)
        resettime_clk= self.qb_pars['resettime']['qubit']
        n_avg = self.qb_pars['n_avg']
        amp_r_scale = self.qb_pars['amp_r_scale']
        # amp_ffl_scale = self.qb_pars['amp_ffl_scale']
        rr_IF = self.qb_pars['rr_IF']
        # ffl_IF = self.qb_pars['ffl_IF']
        with program() as prog:
            update_frequency('rr', rr_IF) 
            # update_frequency('ffl', ffl_IF) 
            n, a, I, Q = declare_vars([int, fixed, fixed, fixed])
        
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            
            with for_(n, 0, n < n_avg, n + 1):
                save(n, n_stream)
                with for_(*from_array(t, t_arr)):
                    play("X90", "qubit")
                    wait(t, "qubit")
                    play("X90", "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None, *res_demod(I, Q))
                    wait(resettime_clk, "qubit")
                    save(I, I_stream)
                    save(Q, Q_stream)
        
            with stream_processing():
                I_stream.buffer(len(t_arr)).average().save("I")
                Q_stream.buffer(len(t_arr)).average().save("Q")
                n_stream.save('n')   

        return prog
    
    def wigner_tomography(self,):
        n_avg = self.qb_pars['n_avg']
        n_points = self.seq_pars['n_points']
        cavity_element = self.qb_pars['cavity_element']
        revival_time = self.qb_pars['revival_time']
        reset_time = self.qb_pars['resettime'][cavity_element]
        alpha = np.linspace(-2, 2, n_points)
        amp_dis = list(-alpha / np.sqrt(2 * np.pi) / 4)
        # ground = declare(int)
        # excited = declare(int)
        with program() as prog:
            n, r, i, I, Q = declare_vars([int, int, int, fixed, fixed, fixed])
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            amp_dis = declare(fixed, value=amp_displace)
            ground_st = declare_stream()
            excited_st = declare_stream()

            with for_(r, 0, r < n_points, r + 1):
                with for_(i, 0, i < n_points, i + 1):
                    assign(ground, 0)
                    assign(excited, 0)
                    with for_(n, 0, n < n_avg, n + 1):
                        # Displace the cavity
                        play("displace" * amp(amp_dis[r], 0, 0, amp_dis[i]), cavity_element)
                        align(cavity_element, "qubit")
                        # The Ramsey sequence with idle time set to pi / chi
                        play("x90", "qubit")
                        wait(revival_time, "qubit")
                        play("x90", "qubit")
                        # Readout the resonator
                        align("qubit", "resonator")
                        measure(
                            "readout",
                            "resonator",
                            None,
                            dual_demod.full("cos", "out1", "sin", "out2", I),
                            dual_demod.full("minus_sin", "out1", "cos", "out2", Q),
                        )
                        # Single shot detection and ground/excited state assignment
                        # with if_(I < threshold):
                        #     assign(ground, ground + 1)
                        # with else_():
                        #     assign(excited, excited + 1)
                        # wait and let all elements relax
                        wait(reset_time, cavity_element, "qubit", "resonator")
                    save(I, I_st)
                    save(Q, Q_st)

        with stream_processing():
            ground_st.buffer(n_points, n_points).save("ground")
            excited_st.buffer(n_points, n_points).save("excited")


        # def make_qubit_reset_sequence():
        #     ##start with a scrambling pulse of length saturation_dur and amplitude scrambling amp on the qubit and then play (or dont play) ffl pulse for length ffl_len and then do rabi.
        #     resettime_clk= clk(qb.pars['qubit_resettime'])
        #     n_avg = self.exp_dict['n_avg']
        #     tmin = clk(tmin)
        #     tmax = clk(tmax)
        #     dt = clk(dt)
        #     amp_ffl_scale = self.exp_dict['amp_ffl_scale']
        #     amp_q_scaling=1.
        #     var_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
        #     with program() as prog:
        #         # n = declare(int)
        #         # I=declare(fixed)
        #         # I_st = declare_stream()
        #         # Q = declare(fixed)
        #         # Q_st = declare_stream()
        #         # n_st=declare_stream()
        #         n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])

        #         I_st, Q_st, n_st = qb.declare_streams(stream_num=3)

        #         update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO'])) # sets the IF frequency of the qubit
        #         update_frequency('rr', qb.pars['rr_IF']) 
        #         update_frequency('ffl', qb.pars['ffl_IF']) 

        #         with for_(n, 0, n < n_avg, n + 1):

        #             save(n, n_st)
        #             with for_each_(t,var_arr):
        #                 play("const" * amp(scrambling_amp), "qubit", duration = clk(saturation_dur)) #play scrambling pulse
        #                 align("qubit", "ffl")
        #                 play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
        #                 align("qubit","ffl")
        #                 play("pi" * amp(amp_q_scaling), "qubit", duration=t)
        #                 align()
        #                 measure("readout", "rr", None, *qb.res_demod(I, Q)) 
        #                 save(I, I_st)
        #                 save(Q, Q_st)
        #                 wait(resettime_clk,"qubit")
                            
        #         with stream_processing():
        #             I_st.buffer(len(var_arr)).average().save("I")
        #             Q_st.buffer(len(var_arr)).average().save("Q")
        #             n_st.save('n')
                    
        
        # elif self.name == 'cavity-reset':
        #     tmin = clk(tmin)
        #     tmax = clk(tmax)
        #     dt = clk(dt)
        #     t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
        #     resettime_clk= clk(qb.pars['qubit_resettime'])
        #     n_avg = self.exp_dict['n_avg']
        #     amp_r_scale = self.exp_dict['amp_r_scale']
        #     amp_ffl_scale = self.exp_dict['amp_ffl_scale']
        #     rr_IF = qb.pars['rr_IF']
        #     ffl_IF = qb.pars['ffl_IF']
        #     qb_IF=qb.pars['qubit_IF']+detuning
        #     with program() as prog:
        #         update_frequency('rr', rr_IF) 
        #         update_frequency('ffl', ffl_IF) 
        #         update_frequency('qubit' ,qb_IF)
        #         n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
        #         I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
        #         with for_(n, 0, n < n_avg, n + 1):
        #             save(n, n_stream)
        #             with for_each_(t, t_arr):
        #                 play("readout"*amp(amp_r_scale), "rr")
        #                 align('ffl','rr')
        #                 #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
        #                 play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
        #                 wait(clk(80))
        #                 align('ffl','qubit')
        #                 play("pi_half", "qubit")
        #                 with if_(t>0):
        #                 #     #align("ffl","qubit")
        #                     wait(t)
        #                 #     #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=t)
        #                 # #     
        #                 # #     #wait(t+clk(20))
        #                 play("pi", "qubit")
        #                 with if_(t>0):
        #                     wait(t, "qubit")
                        
        #                 play("pi_half", "qubit")
        #                 align("qubit","rr")
        #                 # play("pi","qubit")
        #                 # align("qubit","rr")
        #                 # with if_(t>0):
        #                 #     play('gaussian'*amp(amp_ffl_scale), "ffl", duration=t)
        #                 # align("ffl","rr")
        #                 measure("readout", "rr", None, *qb.res_demod(I, Q))
        #                 save(I, I_stream)
        #                 save(Q, Q_stream)
        #                 wait(resettime_clk, 'rr')
        #         with stream_processing():
        #             I_stream.buffer(len(t_arr)).average().save("I")
        #             Q_stream.buffer(len(t_arr)).average().save("Q")
        #             n_stream.save('n')
        
        
        # elif self.name == 'cavity-cooling':
        #     tmin = clk(tmin)
        #     tmax = clk(tmax)
        #     dt = clk(dt)
        #     t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
        #     resettime_clk= clk(qb.pars['qubit_resettime'])
        #     n_avg = self.exp_dict['n_avg']
        #     amp_r_scale = self.exp_dict['amp_r_scale']
        #     amp_ffl_scale = self.exp_dict['amp_ffl_scale']
        #     rr_IF = qb.pars['rr_IF']
        #     ffl_IF = qb.pars['ffl_IF']
        #     qb_IF=qb.pars['qubit_IF']+detuning
        #     with program() as prog:
        #         update_frequency('rr', rr_IF) 
        #         update_frequency('ffl', ffl_IF) 
        #         update_frequency('qubit' ,qb_IF)
        #         n, t, I, Q= qb.declare_vars([int, int, fixed, fixed])
        #         I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
        #         with for_(n, 0, n < n_avg, n + 1):
        #             save(n, n_stream)
        #             with for_each_(t, t_arr):
        #                 #play("readout"*amp(amp_r_scale), "rr")
        #                 #align('ffl','rr')
        #                 #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
        #                 #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len), condition= with_ffl==True)
        #                 #wait(clk(200))
        #                 #align('ffl','qubit')
        #                 #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=2*t+)
        #                 play('const'*amp(amp_ffl_scale), "ffl", duration=2*t+100)
        #                 wait(40,"qubit")
        #                 play("pi_half", "qubit")
        #                 align('rr','qubit')
        #                 play("readout"*amp(amp_r_scale), "rr", duration=2*t+20)
        #                 with if_(t>0):
        #                     wait(t,"qubit")
        #                 play("pi", "qubit")
        #                 with if_(t>0):
        #                     wait(t,"qubit")
        #                 play("pi_half", "qubit")
        #                 #align("qubit","rr")
        #                 align("ffl","rr")
        #                 measure("readout", "rr", None, *qb.res_demod(I, Q))
        #                 save(I, I_stream)
        #                 save(Q, Q_stream)
        #                 wait(resettime_clk)
        #         with stream_processing():
        #             I_stream.buffer(len(t_arr)).average().save("I")
        #             Q_stream.buffer(len(t_arr)).average().save("Q")
        #             n_stream.save('n')
                    
                    
        # elif self.name == 'cavity-cooling-ramsey':
        #     tmin = clk(tmin)
        #     tmax = clk(tmax)
        #     dt = clk(dt)
        #     t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
        #     resettime_clk= clk(qb.pars['qubit_resettime'])
        #     n_avg = self.exp_dict['n_avg']
        #     amp_r_scale = self.exp_dict['amp_r_scale']
        #     amp_ffl_scale = self.exp_dict['amp_ffl_scale']
        #     rr_IF = qb.pars['rr_IF']
        #     ffl_IF = qb.pars['ffl_IF']
        #     qb_IF=qb.pars['qubit_IF']+detuning
        #     with program() as prog:
        #         update_frequency('rr', rr_IF) 
        #         update_frequency('ffl', ffl_IF) 
        #         update_frequency('qubit' ,qb_IF+detuning)
        #         n, t, I, Q= qb.declare_vars([int, int, fixed, fixed])
        #         I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
        #         with for_(n, 0, n < n_avg, n + 1):
        #             save(n, n_stream)
        #             with for_each_(t, t_arr):
        #                 #play("readout"*amp(amp_r_scale), "rr")
        #                 #align('ffl','rr')
        #                 #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
        #                 #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len), condition= with_ffl==True)
        #                 #wait(clk(200))
        #                 #align('ffl','qubit')
        #                 #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=2*t+)
        #                 play('const'*amp(amp_ffl_scale), "ffl", duration=t+90)
        #                 wait(40,"qubit")
        #                 play("pi_half", "qubit")
        #                 align('rr','qubit')
        #                 play("readout"*amp(amp_r_scale), "rr", duration=t)
        #                 with if_(t>0):
        #                     wait(t,"qubit")
        #                 play("pi_half", "qubit")
        #                 #align("qubit","rr")
        #                 align("ffl","rr")
        #                 measure("readout", "rr", None, *qb.res_demod(I, Q))
        #                 save(I, I_stream)
        #                 save(Q, Q_stream)
        #                 wait(resettime_clk)
        #         with stream_processing():
        #             I_stream.buffer(len(t_arr)).average().save("I")
        #             Q_stream.buffer(len(t_arr)).average().save("Q")
        #             n_stream.save('n')
        

    
                    
        # return prog
        
    
    
    def simulate_sequence(self,qb, duration):
        duration = clk(duration)
        qmm = QuantumMachinesManager(host=host, port=port)
        prog = self.make_sequence(qb)
        job = qmm.simulate(qb.config, prog, SimulationConfig(duration=duration))
        samples = job.get_simulated_samples()
        samples.con1.plot()
        qmm.close_all_quantum_machines()
        return samples

    
# def main():
    # qb = dissipator('diss08_11a',device_name='diss08_11a')
    # qb.update_value('ffl_freq', 3.07e9)
    # qb.update_value('ffl_IF', 350e6)
    # qb.update_value('ffl_LO', qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
    #qb.update_value('rr_pulse_len_in_clk', 20)
    # seq = sequence('qb-reset', IF_min=45e6,IF_max=54e6,df=0.1e6, res_ringdown_time=int(40))
    # samples = seq.simulate_sequence(qb, duration=3000)
    #qb.update_value('rr_pulse_len_in_clk', 500)
    
if __name__ == "__main__":
    main()

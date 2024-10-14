# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:00:18 2023

@author: lfl
"""
from qm.qua import *
from qualang_tools.loops import from_array
from qm import LoopbackInterface, QuantumMachinesManager, SimulationConfig
import numpy as np
from Utilities import clk
from helper_functions import res_demod, declare_vars, declare_streams

class sequence:
    def __init__(self, qb, name, **kwargs):
        self.name = name
        self.qb_pars = qb.pars
        self.seq_pars = {key: kwargs.get(key) for key in kwargs.keys()}

    def single_tone_spectroscopy(self):
        IF_min, IF_max, df = self.seq_pars['IF_min'], self.seq_pars['IF_max'], self.seq_pars['df']
        freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
        target_elem = self.seq_pars['target_elem']
        IQ_rotation = self.qb_pars['IQ_rotation']
        phase_prefactor = self.qb_pars['tof'] * 1e-3 * 0.000064
        initial_frame_rotation = self.qb_pars['tof'] * IF_min / 1e9

        with program() as prog:
            n, I, Q, f = declare_vars([int, fixed, fixed, int])
            I_st, Q_st, n_stream = declare_streams(stream_num=3)
            phase_shift = declare(fixed)
            df = declare(int, value=int(df))

            with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                save(n, n_stream)
                reset_frame('rr')
                frame_rotation_2pi(IQ_rotation + initial_frame_rotation, 'rr')
                with for_(*from_array(f, freqs)):
                    update_frequency("rr", f, keep_phase=False)
                    assign(phase_shift, Cast.mul_fixed_by_int(phase_prefactor, (df >> 6)))
                    frame_rotation_2pi(phase_shift, 'rr')
                    wait(self.qb_pars['resettime'][target_elem], target_elem)
                    measure("readout", "rr", None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                    save(I, I_st)
                    save(Q, Q_st)

            with stream_processing():
                I_st.buffer(len(freqs)).average().save('I')
                Q_st.buffer(len(freqs)).average().save('Q')
                n_stream.save('n')

        return prog

    def two_tone_spectroscopy(self):
        freqs = np.arange(self.seq_pars['IF_min'], self.seq_pars['IF_max'] + self.seq_pars['df']/2, self.seq_pars['df'], dtype=int)
        target_elem, readout_elem = self.seq_pars['target_res'], self.seq_pars['readout_res']
        

        with program() as prog:
            n, I, Q, f = declare_vars([int, fixed, fixed, int])
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            if self.seq_pars['on_off']:
                I_b, Q_b, I_tot, Q_tot = declare_vars([fixed, fixed, fixed, fixed])

            with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                save(n, n_stream)
                with for_(*from_array(f, freqs)):
                    update_frequency(target_elem, f)
                    if self.seq_pars['on_off']:
                        measure("readout", readout_elem, None, *res_demod(I_b, Q_b, switch_weights=self.qb_pars['switch_weights']))
                        align(readout_elem, target_elem)
                        wait(self.qb_pars['resettime'][readout_elem], readout_elem)
                        align(readout_elem, target_elem)

                    play('const' * amp(self.seq_pars['amp_q_scaling']), target_elem, duration=self.seq_pars['saturation_dur'])
                    align(target_elem, readout_elem)
                    measure("readout", readout_elem, None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                    wait(self.qb_pars['resettime'][target_elem], target_elem)

                    if self.seq_pars['on_off']:
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

    def two_tone_spectroscopy_amp_sweep(self):
        freqs = np.arange(self.seq_pars['IF_min'], self.seq_pars['IF_max'] + self.seq_pars['df']/2, self.seq_pars['df'], dtype=int)
        target_elem, readout_elem = self.seq_pars['target_res'], self.seq_pars['readout_res']
        amp_arr = np.arange(self.seq_pars['amin'], self.seq_pars['amax'] + self.seq_pars['da']/2, self.seq_pars['da'], dtype=float)

        with program() as prog:
            n, I, Q = declare_vars([int, fixed, fixed])
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            if self.seq_pars['on_off']:
                I_b, Q_b, I_tot, Q_tot = declare_vars([fixed, fixed, fixed, fixed])
            f, a = declare_vars([int, fixed])

            with for_(*from_array(a, amp_arr)):
                with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                    save(n, n_stream)
                    with for_(*from_array(f, freqs)):
                        update_frequency(target_elem, f)
                        if self.seq_pars['on_off']:
                            measure("readout", readout_elem, None, *res_demod(I_b, Q_b, switch_weights=self.qb_pars['switch_weights']))
                            wait(self.qb_pars['resettime'][readout_elem], readout_elem)
                            align(readout_elem, target_elem)

                        play('const' * amp(a), target_elem, duration=self.seq_pars['saturation_dur'])
                        align(target_elem, readout_elem)
                        measure("readout", readout_elem, None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                        wait(self.qb_pars['resettime'][target_elem], target_elem)

                        if self.seq_pars['on_off']:
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

    def power_rabi(self):
        amin, amax, da = self.seq_pars['amin'], self.seq_pars['amax'], self.seq_pars['da']
        amps = np.arange(amin, amax + da/2, da, dtype=float)

        with program() as prog:
            n, a, I, Q = declare_vars([int, fixed, fixed, fixed])
            I_stream, Q_stream, n_stream, a_stream = declare_streams(stream_num=4)

            with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                save(n, n_stream)
                with for_(*from_array(a, amps)):
                    play("gauss" * amp(a), "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                    wait(self.qb_pars['resettime']['qubit'], "qubit")
                    save(I, I_stream)
                    save(Q, Q_stream)
                    save(a, a_stream)

            with stream_processing():
                I_stream.buffer(len(amps)).average().save("I")
                Q_stream.buffer(len(amps)).average().save("Q")
                a_stream.buffer(len(amps)).save('amps')
                n_stream.save('n')

        return prog

    def time_rabi(self):
        tmin, tmax, dt = self.seq_pars['tmin'], self.seq_pars['tmax'], self.seq_pars['dt']
        t_arr = np.arange(tmin, tmax + dt/2, dt, dtype=int)

        with program() as prog:
            n, t, I, Q = declare_vars([int, int, fixed, fixed])
            I_stream, Q_stream, n_stream, t_stream = declare_streams(stream_num=4)

            with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                save(n, n_stream)
                with for_(*from_array(t, t_arr)):
                    play("gauss", "qubit",duration=t)
                    align("qubit", "rr")
                    measure("readout", "rr", None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                    wait(self.qb_pars['resettime']['qubit'], "qubit")
                    save(I, I_stream)
                    save(Q, Q_stream)
                    save(t, t_stream)

            with stream_processing():
                I_stream.buffer(len(t_arr)).average().save("I")
                Q_stream.buffer(len(t_arr)).average().save("Q")
                t_stream.buffer(len(t_arr)).save('times')
                n_stream.save('n')

        return prog
    
    def ramsey(self):
        tmin, tmax, dt = self.seq_pars['tmin'], self.seq_pars['tmax'], self.seq_pars['dt']
        t_arr = np.arange(tmin, tmax + dt/2, dt, dtype=int)

        with program() as prog:
            n, t, I, Q = declare_vars([int, int, fixed, fixed])
            I_stream, Q_stream, n_stream, t_stream = declare_streams(stream_num=4)

            with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                save(n, n_stream)
                with for_(*from_array(t, t_arr)):
                    play("X90", "qubit")
                    wait(t, "qubit")
                    play("X90", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                    wait(self.qb_pars['resettime']['qubit'], "qubit")
                    save(I, I_stream)
                    save(Q, Q_stream)
                    save(t, t_stream)

            with stream_processing():
                I_stream.buffer(len(t_arr)).average().save("I")
                Q_stream.buffer(len(t_arr)).average().save("Q")
                t_stream.buffer(len(t_arr)).save('times')
                n_stream.save('n')

        return prog

    def T1(self):
        tmin, tmax, dt = self.seq_pars['tmin'], self.seq_pars['tmax'], self.seq_pars['dt']
        t_arr = np.arange(tmin, tmax + dt/2, dt, dtype=int)

        with program() as prog:
            n, t, I, Q = declare_vars([int, int, fixed, fixed])
            I_stream, Q_stream, n_stream, t_stream = declare_streams(stream_num=4)

            with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                save(n, n_stream)
                with for_(*from_array(t, t_arr)):
                    play("X180", "qubit")
                    wait(t, "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                    wait(self.qb_pars['resettime']['qubit'], "qubit")
                    save(I, I_stream)
                    save(Q, Q_stream)
                    save(t, t_stream)

            with stream_processing():
                I_stream.buffer(len(t_arr)).average().save("I")
                Q_stream.buffer(len(t_arr)).average().save("Q")
                t_stream.buffer(len(t_arr)).save('times')
                n_stream.save('n')

        return prog

    def echo(self):
        tmin, tmax, dt = self.seq_pars['tmin'], self.seq_pars['tmax'], self.seq_pars['dt']
        t_arr = np.arange(tmin, tmax + dt/2, dt, dtype=int)

        with program() as prog:
            n, t, I, Q = declare_vars([int, int, fixed, fixed])
            I_stream, Q_stream, n_stream, t_stream = declare_streams(stream_num=4)

            with for_(n, 0, n < self.qb_pars['n_avg'], n + 1):
                save(n, n_stream)
                with for_(*from_array(t, t_arr)):
                    play("X90", "qubit")
                    wait(t, "qubit")
                    play("X180", "qubit")
                    wait(t, "qubit")
                    play("X90", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *res_demod(I, Q, switch_weights=self.qb_pars['switch_weights']))
                    wait(self.qb_pars['resettime']['qubit'], "qubit")
                    save(I, I_stream)
                    save(Q, Q_stream)
                    save(t, t_stream)

            with stream_processing():
                I_stream.buffer(len(t_arr)).average().save("I")
                Q_stream.buffer(len(t_arr)).average().save("Q")
                t_stream.buffer(len(t_arr)).save('times')
                n_stream.save('n')

        return prog

    def wigner_tomography(self):
        n_avg = self.qb_pars['n_avg']
        n_points = self.seq_pars['n_points']
        cavity_element = self.qb_pars['cavity_element']
        revival_time = self.qb_pars['revival_time']
        reset_time = self.qb_pars['resettime'][cavity_element]
        alpha = np.linspace(-2, 2, n_points)
        amp_dis = list(-alpha / np.sqrt(2 * np.pi) / 4)

        with program() as prog:
            n, r, i, I, Q = declare_vars([int, int, int, fixed, fixed])
            I_stream, Q_stream, n_stream = declare_streams(stream_num=3)
            amp_dis = declare(fixed, value=amp_displace)
            ground_st = declare_stream()
            excited_st = declare_stream()

            with for_(r, 0, r < n_points, r + 1):
                with for_(i, 0, i < n_points, i + 1):
                    assign(ground, 0)
                    assign(excited, 0)
                    with for_(n, 0, n < n_avg, n + 1):
                        play("displace" * amp(amp_dis[r], 0, 0, amp_dis[i]), cavity_element)
                        align(cavity_element, "qubit")
                        play("x90", "qubit")
                        wait(revival_time, "qubit")
                        play("x90", "qubit")
                        align("qubit", "resonator")
                        measure("readout", "resonator", None, dual_demod.full("cos", "out1", "sin", "out2", I), dual_demod.full("minus_sin", "out1", "cos", "out2", Q))
                        wait(reset_time, cavity_element, "qubit", "resonator")
                    save(I, I_st)
                    save(Q, Q_st)

        with stream_processing():
            ground_st.buffer(n_points, n_points).save("ground")
            excited_st.buffer(n_points, n_points).save("excited")

    def simulate_sequence(self, qb, duration):
        duration = clk(duration)
        qmm = QuantumMachinesManager(host=host, port=port)
        prog = self.make_sequence(qb)
        job = qmm.simulate(qb.config, prog, SimulationConfig(duration=duration))
        samples = job.get_simulated_samples()
        samples.con1.plot()
        qmm.close_all_quantum_machines()
        return samples

if __name__ == "__main__":
    main()

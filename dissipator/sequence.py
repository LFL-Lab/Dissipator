# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:00:18 2023

@author: lfl
"""
from qm import generate_qua_script
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm import SimulationConfig
import instrument_init as inst
from dissipator import *
from Utilities import clk
host='10.71.0.57'
port='9510'
class sequence():
    
    def __init__(self, name, n_avg=100, amp_r_scale=1, amp_ffl_scale=1, **kwargs):
        self.name = name
        self.exp_dict = {'n_avg': n_avg,
                         'amp_r_scale': amp_r_scale,
                         'amp_ffl_scale': amp_ffl_scale}
        for key in kwargs.keys():
            self.exp_dict[key] = kwargs.get(key)
        
        
    def make_sequence(self, qb, tmin = 16, tmax=1e3, dt=8, scrambling_amp=1, ffl_len=0.5e3, saturation_dur=1e3, with_ffl='True',):
        if self.name == 'ringdown_drive_on':
            tmin = clk(tmin)
            tmax = clk(tmax)
            dt = clk(dt)
            t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
            resettime_clk= clk(qb.pars['rr_resettime'])
            n_avg = self.exp_dict['n_avg']
            amp_r_scale = self.exp_dict['amp_r_scale']
            amp_ffl_scale = self.exp_dict['amp_ffl_scale']
            rr_IF = qb.pars['rr_IF']
            ffl_IF = qb.pars['ffl_IF']
            with program() as prog:
                update_frequency('rr', rr_IF) 
                update_frequency('ffl', ffl_IF) 
                n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
            
                I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
                
                with for_(n, 0, n < n_avg, n + 1):
                    save(n, n_stream)
                    with for_each_(t, t_arr):
                        with if_(t==0):
                            play("readout"*amp(amp_r_scale), "rr")
                            measure("void", "rr", None,*qb.res_demod(I,Q))
                            wait(resettime_clk, "rr")
                            save(I, I_stream)
                            save(Q, Q_stream)
                        with else_():
                            play("readout"*amp(amp_r_scale), "rr")
                            align("ffl", "rr")
                            play('const'*amp(amp_ffl_scale), "ffl", duration=t)
                            align("ffl", "rr")
                            measure("void", "rr", None,*qb.res_demod(I,Q))
                            wait(resettime_clk, "rr")
                            save(I, I_stream)
                            save(Q, Q_stream)
            
                with stream_processing():
                    I_stream.buffer(len(t_arr)).average().save("I")
                    Q_stream.buffer(len(t_arr)).average().save("Q")
                    n_stream.save('n')

        elif self.name == 'rr_spec':
            n_avg = self.exp_dict['n_avg']
            IF_min = self.exp_dict['IF_min']
            IF_max = self.exp_dict['IF_max']
            df = self.exp_dict['df']
            freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
            rr_pulse_len_in_clk = qb.pars['rr_pulse_len_in_clk']
            res_ringdown_time = clk(qb.pars['rr_resettime'])
            ### QUA code ###
            with program() as prog:
             
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
                
                        align('rr', 'ffl')
                
                        measure("readout", "rr", None,*qb.res_demod(I, Q))
                        save(I, I_st)
                        save(Q, Q_st)

                    save(n,n_stream)

                with stream_processing():
                    I_st.buffer(len(freqs)).average().save('I')
                    Q_st.buffer(len(freqs)).average().save('Q')
                    n_stream.save('n')
					
		elif self.name == 'qb-reset-with-scramb':
            ##start with a scrambling pulse of length saturation_dur and amplitude scrambling amp on the qubit and then play (or dont play) ffl pulse for length ffl_len and then do rabi.
            resettime_clk= clk(qb.pars['qubit_resettime'])
			n_avg = self.exp_dict['n_avg']
			tmin = clk(tmin)
            tmax = clk(tmax)
            dt = clk(dt)
			amp_ffl_scale = self.exp_dict['amp_ffl_scale']
			var_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
            with program() as prog:
				n = declare(int)
                I=declare(fixed)
                I_st = declare_stream()
               	Q = declare(fixed)
               	Q_st = declare_stream()
                n_st=declare_stream()
                #n, t, I, Q = self.declare_vars([int, int, fixed, fixed])

                #I_stream, Q_stream, n_stream = self.declare_streams(stream_num=3)

                update_frequency('qubit', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning) # sets the IF frequency of the qubit
				update_frequency('rr', rr_IF) 
                update_frequency('ffl', ffl_IF) 

                with for_(n, 0, n < n_avg, n + 1):
					
                    save(n, n_st)
                    with for_each_(t,var_arr):
                        with if_(with_ffl == 'True'):
                            with if_(t==0):
                                play("const" * amp(scrambling_amp), "qubit", duration = clk(saturation_dur)) #play scrambling pulse
    							align("ffl", "qubit")
                                play('const'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
    							align("rr","ffl")
                                measure("readout", "rr", None, *self.res_demod(I, Q))                     
                                # save(t,t_stream)
                                save(I, I_st)
                                save(Q, Q_st)
                                wait(resettime_clk,"qubit")
                            with else_():
    							play("const" * amp(scrambling_amp), "qubit", duration = clk(saturation_dur)) #play scrambling pulse
    							align("ffl", "qubit")
                                play('const'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
    							align("qubit","ffl")
                                play("pi" * amp(amp_q_scaling), "qubit", duration=t)
                                align("qubit", "rr")
                                measure("readout", "rr", None, *self.res_demod(I, Q))
                                # save(t,t_stream)
                                save(I, I_st)
                                save(Q, Q_st)
                                wait(resettime_clk,"qubit")
                        with else_():
                            with if_(t==0):
                                play("const" * amp(scrambling_amp), "qubit", duration = clk(saturation_dur)) #play scrambling pulse
    							align("rr", "qubit")
                                measure("readout", "rr", None, *self.res_demod(I, Q))                     
                                # save(t,t_stream)
                                save(I, I_st)
                                save(Q, Q_st)
                                wait(resettime_clk,"qubit")
                            with else_()
                                play("const" * amp(scrambling_amp), "qubit", duration = clk(saturation_dur)) #play scrambling pulse
                                play("pi" * amp(amp_q_scaling), "qubit", duration=t)
                                align("qubit", "rr")
                                measure("readout", "rr", None, *self.res_demod(I, Q))
                                # save(t,t_stream)
                                save(I, I_st)
                                save(Q, Q_st)
                                wait(resettime_clk,"qubit")
                            

                with stream_processing():
                    I_st.buffer(len(var_arr)).average().save("I")
                    Q_st.buffer(len(var_arr)).average().save("Q")
                    n_st.save('n')
            
            elif self.name == 'qb-reset-without-scramb':
                #start with ffl pulse for length ffl_len, and then do rabi. includes a version where no ffl pulse is present (equivalent to pure rabi.)
                resettime_clk= clk(qb.pars['qubit_resettime'])
    			n_avg = self.exp_dict['n_avg']
    			tmin = clk(tmin)
                tmax = clk(tmax)
                dt = clk(dt)
    			amp_ffl_scale = self.exp_dict['amp_ffl_scale']
    			var_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
                with program() as prog:
    				n = declare(int)
                    I=declare(fixed)
                    I_st = declare_stream()
                   	Q = declare(fixed)
                   	Q_st = declare_stream()
                    n_st=declare_stream()
                    #n, t, I, Q = self.declare_vars([int, int, fixed, fixed])

                    #I_stream, Q_stream, n_stream = self.declare_streams(stream_num=3)

                    update_frequency('qubit', (self.pars['qubit_freq']-self.pars['qubit_LO']) + detuning) # sets the IF frequency of the qubit
    				update_frequency('rr', rr_IF) 
                    update_frequency('ffl', ffl_IF) 

                    with for_(n, 0, n < n_avg, n + 1):
    					
                        save(n, n_st)
                        with for_each_(t,var_arr):
                            with if_(with_ffl == 'True'):
                                with if_(t==0):
                                    play('const'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
        							align("rr","ffl")
                                    measure("readout", "rr", None, *self.res_demod(I, Q))                     
                                    # save(t,t_stream)
                                    save(I, I_st)
                                    save(Q, Q_st)
                                    wait(resettime_clk,"qubit")
                                with else_():
                                    play('const'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
        							align("qubit","ffl")
                                    play("pi" * amp(amp_q_scaling), "qubit", duration=t)
                                    align("qubit", "rr")
                                    measure("readout", "rr", None, *self.res_demod(I, Q))
                                    # save(t,t_stream)
                                    save(I, I_st)
                                    save(Q, Q_st)
                                    wait(resettime_clk,"qubit")
                            with else_():
                                with if_(t==0):
                                    measure("readout", "rr", None, *self.res_demod(I, Q))                     
                                    # save(t,t_stream)
                                    save(I, I_st)
                                    save(Q, Q_st)
                                    wait(resettime_clk,"qubit")
                                with else_()
                                    play("pi" * amp(amp_q_scaling), "qubit", duration=t)
                                    align("qubit", "rr")
                                    measure("readout", "rr", None, *self.res_demod(I, Q))
                                    # save(t,t_stream)
                                    save(I, I_st)
                                    save(Q, Q_st)
                                    wait(resettime_clk,"qubit")
                                

                    with stream_processing():
                        I_st.buffer(len(var_arr)).average().save("I")
                        Q_st.buffer(len(var_arr)).average().save("Q")
                        n_st.save('n')                 
                    
        return prog
        
    
    
    def simulate_sequence(self,qb, duration):
        duration = clk(duration)
        qmm = QuantumMachinesManager(host=host, port=port)
        prog = self.make_sequence(qb)
        job = qmm.simulate(qb.config, prog, SimulationConfig(duration=duration))
        samples = job.get_simulated_samples()
        samples.con1.plot(analog_ports=['1','2','3', '4', '5', '6'])
        qmm.close_all_quantum_machines()
        return samples
    
def main():
    device = 'diss09_6024'
    qb = dissipator(device, device_name=device)
    # qb.update_value('ffl_freq', 3.07e9)
    # qb.update_value('ffl_IF', 350e6)
    # qb.update_value('ffl_LO', qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
    qb.update_value('rr_pulse_len_in_clk', 20)
    seq = sequence('rr_spec', IF_min=45e6,IF_max=54e6,df=0.1e6, res_ringdown_time=int(40))
    samples = seq.simulate_sequence(qb, duration=4000)
    qb.update_value('rr_pulse_len_in_clk', 500)
    
if __name__ == "__main__":
    main()

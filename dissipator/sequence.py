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
host='10.71.0.57'
port='9510'
class sequence():
    
    def __init__(self, name, n_avg=100, amp_r_scale=1, amp_ffl_scale=1):
        self.name = name
        self.exp_dict = {'n_avg': n_avg,
                         'amp_r_scale': amp_r_scale,
                         'amp_ffl_scale': amp_ffl_scale}
        
        
        
    def make_sequence(self, qb, tmin = 16, tmax=1e3, dt=8):
        if self.name == 'ringdown_drive_on':
            tmin = clk(tmin)
            tmax = clk(tmax)
            dt = clk(dt)
            t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
            resettime_clk= 100
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
                    with for_(*from_array(t,t_arr)):
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
                
        return prog
    
    def simulate_sequence(self,qb, duration):
        duration = clk(duration)
        qmm = QuantumMachinesManager(host=host, port=port)
        prog = self.make_sequence(qb)
        job = qmm.simulate(qb.config, prog, SimulationConfig(duration=duration))
        samples = job.get_simulated_samples()
        samples.con2.plot(analog_ports=['1','2','3', '4', '5', '6'])
        qmm.close_all_quantum_machines()
        return samples
        
qb = dissipator('diss09', device_name='diss09_5578')
qb.update_value('ffl_freq', 3.07e9)
qb.update_value('ffl_IF', 350e6)
qb.update_value('ffl_LO', qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
seq = sequence('ringdown_drive_on')
samples = seq.simulate_sequence(qb, duration=4000)

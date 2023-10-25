# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 15:04:10 2023

@author: lfl
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 16:00:14 2022

@author: lfl
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 15:28:33 2022

@author: lfl
"""

from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
#from configuration_NIST_Q2 import *
import matplotlib.pyplot as plt
import numpy as np

###################
# The QUA program #
###################

f_min = -200e6
f_max = -205e6
df = 0.2e6

freqs = np.arange(f_min, f_max + df/2, df)

with program() as IQ_blobs:

    n = declare(int)
    I = declare(fixed)
    I_st = declare_stream()
    Q = declare(fixed)
    Q_st = declare_stream()
    I_exc = declare(fixed)
    Q_exc = declare(fixed)
    I_st_exc = declare_stream()
    Q_st_exc = declare_stream()
    f = declare(int)
    distance = declare(fixed)
    distance_st = declare_stream()
    resettime_clk= clk(self.pars['qubit_resettime'])
    with for_(n, 0, n < 2000, n + 1):
        with for_(f, f_min, f <= f_max, f + df):
            update_frequency('rr', f)
            
            wait(resettime_clk, "qubit")
            align("qubit", "rr")
            measure("readout", "rr", None,
                    dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                    dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
            save(I, I_st)
            save(Q, Q_st)
            
            align('qubit','rr')
            wait(resetTime, "qubit")
            play("pi", "qubit")
            align("qubit", "rr")
            measure("readout", "rr", None,
                    dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I_exc),
                    dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q_exc))
            save(I_exc, I_st_exc)
            save(Q_exc, Q_st_exc)
            
            assign(distance, Math.abs(I_exc - I))
            save(distance, distance_st)
            
    with stream_processing():
        I_st.save_all('I')
        Q_st.save_all('Q')
        I_st_exc.save_all('I_exc')
        Q_st_exc.save_all('Q_exc')
        distance_st.buffer(len(freqs)).average().save('distance')

######################################
# Open Communication with the Server #
######################################
qmm = QuantumMachinesManager()

####################
# Simulate Program #
####################
# simulation_config = SimulationConfig(
#                     duration=5000,
#                     simulation_interface=LoopbackInterface([("con1", 3, "con1", 1), ("con1", 4, "con1", 2)]))

# job = qmm.simulate(config, IQ_blobs, simulation_config)
qm = qmm.open_qm(config)
job = qm.execute(IQ_blobs)
res_handles = job.result_handles
res_handles.wait_for_all_values()
I = res_handles.get("I").fetch_all()['value']
Q = res_handles.get("Q").fetch_all()['value']
I_exc = res_handles.get("I_exc").fetch_all()['value']
Q_exc = res_handles.get("Q_exc").fetch_all()['value']

distance = res_handles.get('distance').fetch_all()

# plt.figure()
# plt.title('IQ_blobs sequence')
# job.get_simulated_samples().con1.plot()
# plt.xlabel("time [ns]")
# plt.ylabel("DAC [V]")

# print(len(I))
# print(len(Q))
# plt.figure()
# plt.plot(I,Q,'.',alpha=0.5)
# plt.axis('equal')
# plt.plot(I_exc,Q_exc,'.',alpha=0.5)
# plt.gca().set_aspect('equal', adjustable='box')

plt.figure()
plt.plot(freqs, distance, '-o')
argmax = int(np.argmax(distance))
freqOpt = freqs[argmax]
plt.plot(freqOpt, distance[argmax], marker='o',color='r')
plt.title(r'$f^*$ = %.3f MHz'%(freqOpt/1e6))
plt.show()# plt.hist2d(I, Q)

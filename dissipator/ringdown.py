# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 18:14:15 2023

@author: lfl
"""

from qubit import *
import instrument_init as inst
import plot_functions as pf
qb = qubit('logical')

tmin = 16
tmax = 2e3
dt = 16
n_avg = 1000
tmin = clk(tmin)
tmax = clk(tmax)
dt = clk(dt)
t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)

inst.set_attenuator(attenuation=qb.pars['rr_atten'])
resettime_clk= clk(qb.pars['qubit_resettime'])

with program() as prog:
    n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])

    I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)

    with for_(n, 0, n < n_avg, n + 1):
        save(n, n_stream)
        with for_(*from_array(t,t_arr)):
            play("readout", "rr")
            wait(t, 'rr')
            measure("void", "rr", None,*qb.res_demod(I,Q))
            wait(resettime_clk, "qubit")
            save(I, I_stream)
            save(Q, Q_stream)

    with stream_processing():
        I_stream.buffer(len(t_arr)).average().save("I")
        Q_stream.buffer(len(t_arr)).average().save("Q")
        n_stream.save('n')

datadict, job = qb.get_results(jobtype = prog, result_names = ["I", "Q"], n_total=n_avg, notify = False)
t_arr = np.array(t_arr)*4/1e3
I = np.array(datadict["I"])
Q = np.array(datadict["Q"])
ydata = np.abs(I+1j*Q)

fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='ringdown',dt=t_arr[-1]*1e-6/len(t_arr))
pf.plot_data(t_arr,ydata,sequence='ringdown',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
             qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
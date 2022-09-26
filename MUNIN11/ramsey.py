# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:52:03 2021

@author: lfl
"""
"""
T2.py: Single qubit T2 decay time
Author: Steven Frankel, Gal Winer - Quantum Machines
Created: 7/11/2020
Created on QUA version: 0.5.138
"""

import matplotlib.pyplot as plt
from config_MUN11 import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import SimulationConfig, LoopbackInterface
import numpy as np
from scipy.optimize import curve_fit
import time
import plot_functions as pf
#from file_utility import *
import pandas as pd

def main(qm,  dformat, bPlot = True, **kwargs):
    t_min = 0
    dt = 16  # timestep
    t_max = 500
    N_max = 5000
    detun = 0
    t_arr = np.arange(t_min, t_max + dt/2 , dt)# Number of timesteps
    resetTime = 5000

    with program() as ramsey:
        I = declare(fixed)# QUA variables declaration
        Q = declare(fixed)
        t = declare(int)  # Sweeping parameter over the set of durations
        Nrep = declare(int)
        phi = declare(fixed,value=-8)
        I_stream = declare_stream()
        Q_stream = declare_stream()
        t_stream = declare_stream()
        Nrep_stream = declare_stream()
        update_frequency('qubit', (qbFreq-qb_LO)+detun)
        # T2
        with for_(Nrep, 0, Nrep < N_max, Nrep + 1):  # Do a 100 times the experiment to obtain statistics
            save(Nrep,Nrep_stream)    
            with for_each_(t, t_arr):
                with if_(t==0):
                    
                    play("pi_half", "qubit")
                    play("pi_half", "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
                    save(t, t_stream)
                    wait(resetTime,"qubit")
                    # assign(phi,phi+0.1)
                
                with else_():
                    play("pi_half", "qubit")
                    wait(t, "qubit")
                    # frame_rotation_2pi(phi, 'qubit')
                    play("pi_half", "qubit")
                    align("qubit","rr")
                    measure("readout", "rr", None,
                            dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                            dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                    save(I, I_stream)
                    save(Q, Q_stream)
                    save(t, t_stream)
                    wait(resetTime,"qubit")
    
        with stream_processing():
            I_stream.buffer(len(t_arr)).average().save("I")
            Q_stream.buffer(len(t_arr)).average().save("Q")
            Nrep_stream.save('n')
            t_stream.save('time')
    
    job = qm.execute(ramsey)
        
    res_handles = job.result_handles
    
    if bPlot:
        I_handle = res_handles.get("I")
        Q_handle = res_handles.get("Q")
        n_handle = res_handles.get('n')
        t_handle = res_handles.get('time')
        I_handle.wait_for_values(1)
        Q_handle.wait_for_values(1)
        n_handle.wait_for_values(1)
        t_handle.wait_for_values(1)
        while(I_handle.is_processing()):
            I = I_handle.fetch_all()
            Q = Q_handle.fetch_all()
            n = n_handle.fetch_all()
            # mag = np.sqrt(I**2 + Q**2)
            # plt.figure()
            # plt.clf()
            plt.plot(t_arr*4, 1e3*I)
            # plt.title('qubit spectroscopy analysis')
            plt.xlabel("time (ns)")
            plt.ylabel("Amplitude (mV)")
            plt.title('n = %d'%(n))
            plt.pause(0.1)
            plt.clf()
    else:
        res_handles.wait_for_all_values()
        I_handle = res_handles.get("I")
        Q_handle = res_handles.get("Q")
        n_handle = res_handles.get('n')
        t_handle = res_handles.get('time')
        
    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    n = n_handle.fetch_all()
    t = t_handle.fetch_all()
        
    if bPlot:
        detuning,T_phi,error = pf.pulse_plot(sequence="ramsey", t_vector = 4*t_arr, y_vector=I,dt=4*dt,qubitDriveFreq=qb_LO + qb_IF +detun,amp_q=pi_amp)

    dataDict = {'time':t_arr * 4,
                 'I': I,
                 'Q': Q}
    
    # if dformat == 'csv':
    #     filename = 'ramsey_detun=%sMHz_%s'%(str(detun/1e6).replace('.', 'd'),datestr)
    #     filename = make_new_filename(filepath + datestr,filename,'csv')
    #     dataDf = pd.DataFrame(dataDict) 
    #     dataDf.to_csv(filepath +datestr +'\\'+ filename)
        
    # elif dformat == 'hdf5':
    #     if 'dname' in kwargs:
    #         dname = kwargs.get('dname')
    #         if 'dindex' in kwargs:
    #             dindex = kwargs.get('dindex')
    #         df = h5py.File(dname, 'a')
    #     else:
    #         expName = 'untitled'
    #         dataname = '%s_tmax=%dcycles_detun=%sMHz_date=%s'%(expName, t_max, detun, datestr)
    #         datapath = filepath + datestr
    #         filename = make_new_filename(datapath, dataname, 'hdf5')
    #         dname = datapath + '\\' + filename
    #         df = h5py.File(dname, 'w')
    #     if 'ramsey/%d' in df.keys():
    #         raise ValueError('ramsey/dindex = %d exists'%(dindex))
    #     else:
    #         gp = df.create_group('ramsey/%d'%(dindex))
    #         for key, value in dataDict.items():
    #             gp.create_dataset(key, data=value)
    #     df.close()
        
        
if __name__ == '__main__':
    qmm = QuantumMachinesManager()  # Reach OPX's IP address
    start = time.time()
    qm = qmm.open_qm(config)
    main(qm, dformat='csv', bPlot=True)
    end = time.time()
    print(end-start)    
    






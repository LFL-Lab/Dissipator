# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:42:53 2022

@author: lfl
"""
sa = sa
from mixer_calibration import get_power,opt_mixer
from qubit import qubit
import instrument_init as inst

ref_H = - 10
ref_L = -45

qb = qubit('qb6')

qb.update_value('qubit_LO', value = 3e9)
qb.update_value('rr_freq', value = 6.70625e9)
qb.update_value('rr_IF', -30e6)
qb.update_value('rr_LO', value = qb.pars['rr_freq'] - qb.pars['rr_IF'])

qb.tof_cal()

'''Qubit mixer calibration'''
qb_lo_leakage = qb.get_power(sa, freq=inst.get_qb_LO(),reference=ref_H,config=True,plot=True)

qb_im_leakage = qb.get_power(sa, freq=inst.get_qb_LO()-50e6,reference=ref_H,config=True,plot=True)

qb_on_power = qb.get_power(sa, freq=inst.get_qb_LO()+50e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
qb.opt_mixer(sa, cal='LO',mode='coarse', freq_span = 1e6, reference=ref_H,element='qubit')
# do a finer sweep in desired range
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine',element='qubit')
# do a coarse sweep to minimize sideband
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference=ref_H, mode='coarse',element='qubit')
# do a finer sweep in desired range
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element='qubit')

'''Readout mixer calibration'''

rr_lo_leakage = qb.get_power(sa, freq=inst.get_rr_LO(),reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

rr_im_leakage = qb.get_power(sa, freq=inst.get_rr_LO()-qb.pars['rr_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

rr_on_power = qb.get_power(sa, freq=inst.get_rr_LO()+qb.pars['rr_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element='rr')
# do a finer sweep in desired range
qb.opt_mixer( sa, cal='LO',  freq_span = 1e6, reference = ref_L, mode='fine',element='rr')
# do a coarse sweep to minimize sideband
qb.opt_mixer( sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element='rr')
# do a finer sweep in desired range
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element='rr')

#%% resonator spectroscopy

I,Q,freqs,job = qb.run_scan(df = 25e3, n_avg = 250, element='resonator', chunksize = 200e6, attenuation=35, lo_min = 4.3e9, lo_max = 4.9e9,
           showprogress=False, res_ringdown_time = int(4e3))

#%% qubit spectroscopy

I,Q,freqs,job = qb.run_scan(df = 100e3, n_avg = 300, element='qubit', chunksize = 200e6, attenuation = 25, lo_min = 4.2e9, lo_max = 4.6e9,
            amp_q_scaling = 0.75, saturation_dur = 20e3, showprogress=False, res_ringdown_time = int(4e3))

#%% Rabi

t_ar, I, Q, job, fitted_pars = qb.pulse_exp(exp='rabi',n_avg = 2000, t0 = 0,tf = 10e3, dt = 500, amp_q_scaling = 1, fit = True, plot = True,
                                     detuning = 0e6, resettime = 400e3)

#%% Ramsey

#%% Echo

# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:42:53 2022

@author: lfl
"""
sa = sa
# from mixer_calibration import get_power,opt_mixer
from qubit import qubit
import instrument_init as inst

ref_H = - 10
ref_L = -45

'''Initialize qubit class'''
qb = qubit('qb3')

'''Update important parameters'''
qb.update_value('qubit_LO', value = 4.48e9)
qb.update_value('qubit_IF', value = 50e6)
qb.update_value('rr_freq', value = 6.70625e9)
qb.update_value('rr_IF', -30e6)
qb.update_value('rr_LO', value = qb.pars['rr_freq'] - qb.pars['rr_IF'])

qb.tof_cal()

'''Qubit mixer calibration
Get leakage power '''
qb.play_pulses()
qb_lo_leakage = qb.get_power(sa, freq=qb.pars['qubit_LO'],reference=ref_H,config=True,plot=True)
qb_im_leakage = qb.get_power(sa, freq=qb.pars['qubit_LO']-qb.pars['qubit_IF'],reference=ref_H,config=True,plot=True)
qb_on_power = qb.get_power(sa, freq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power

'''Optimize mixer'''
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'qubit')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'qubit')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine', element = 'qubit')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'qubit')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'qubit')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element = 'qubit')

'''Readout mixer calibration'''

rr_lo_leakage = qb.get_power(sa, freq=qb.pars['rr_LO'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
rr_im_leakage = qb.get_power(sa, freq=qb.pars['rr_LO']-qb.pars['rr_IF'],span = 1e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
rr_on_power = qb.get_power(sa, freq=qb.pars['rr_LO']+qb.pars['rr_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
qb.opt_mixer( sa, cal='LO',  freq_span = 1e6, reference = ref_L, mode='fine',element='rr')
qb.opt_mixer( sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element='rr')

#%% resonator spectroscopy

I,Q,freqs,job = qb.run_scan(df = 25e3, n_avg = 250, element='resonator', chunksize = 200e6, attenuation=35, lo_min = 4.3e9, lo_max = 4.9e9,
           showprogress=False, res_ringdown_time = int(4e3))

#%% qubit spectroscopy

I,Q,freqs,job = qb.run_scan(df = 50e3, n_avg = 300, element='qubit', chunksize = 200e6, attenuation = 25, lo_min = 4.3e9, lo_max = 4.3e9,
            amp_q_scaling = 0.05, saturation_dur = 20e3, showprogress=False, res_ringdown_time = int(4e3))

#%% Rabi
atten = 20

t_arr, I, Q, job, pi2width = qb.pulse_exp(exp='rabi',n_avg = 200, atten = 22, tmin = 20,tmax = 1e3, dt = 10, amp_q_scaling = 0.1, fit = True, plot = True,
                                     detuning = 0e6, resettime = 400e3)

amps, I, Q, job, pi_amp = qb.power_rabi(a_min = 0.01, a_max = 1, da = 0.005,  check_mix = False, n_avg = 200,numPeriods=2,update_amp=True,
                                            fit = True, plot = True, detuning = 0e6,  resettime = 400e3)

qb.cal_pi_pulse(pi_half_len_target = 20, starting_amp = 0.44, n_avg = 200)


#%% Ramsey
t_arr, I, Q, job, fitted_pars = qb.pulse_exp(exp='ramsey', check_mixers = True, n_avg=200,dt=200,tmin=20,tmax=4e4,detuning=0e3)

#%% Echo

#%% T1

t_arr, I, Q,job,fitted_pars = qb.pulse_exp(exp='T1',n_avg=200,atten=20,tmin=100,tmax=10e4,dt=1000,detuning=0)

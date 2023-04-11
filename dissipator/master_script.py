# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:42:53 2022

@author: lfl
"""
sa = sa
# from mixer_calibration import get_power,opt_mixer

from qubit import qubit
import instrument_init as inst
import plot_functions as pf
import numpy as np
from ringdown_ffl import measure_ringdown_drive_off, measure_ringdown_drive_on
from curvefitgui import curve_fit_gui
from t1_ffl import measure_t1_ffl_off, measure_t1_w_ffl

ref_H = 20
ref_L = -45

'''Initialize qubit class'''
qb = qubit('logical')


'''Update important parameters'''
qb.update_value('ffl_freq', 2.85e9)
qb.update_value('ffl_LO', 2.8e9)
qb.update_value('ffl_IF', 0.05e9)
qb.update_value('qubit_LO', 4.5e9)
qb.update_value('qubit_freq', 4.6383e9)
qb.update_value('qubit_IF', qb.pars['qubit_freq'] - qb.pars['qubit_LO'])

'''Qubit mixer calibration 
Get leakage power '''
qb.play_pulses()
inst.set_qb_LO(qb.pars['qubit_LO'])

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
# set DA to 0 dB attenuation
set_attenuator(0)
get_attenuation()
rr_lo_leakage = qb.get_power(sa, freq=qb.pars['rr_LO'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
rr_im_leakage = qb.get_power(sa, freq=qb.pars['rr_LO']-qb.pars['rr_IF'],span = 1e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
rr_on_power = qb.get_power(sa, freq=qb.pars['rr_LO']+qb.pars['rr_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine',element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element='rr')

'''FFL mixer calibration'''
ref_H = 0
ref_L = -30

qb.play_pulses()
sa = init_sa_by_serial_number(20234492)
inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
qb_lo_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO'],reference=ref_H,config=True,plot=True)
qb_im_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO']-qb.pars['ffl_IF'],reference=ref_H,config=True,plot=True)
qb_on_power = qb.get_power(sa, freq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power

'''Optimize FFL mixer'''
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine', element = 'ffl')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
# qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element = 'ffl')


'''diss mixer calibration'''
inst.set_diss_LO(qb.pars['diss_LO'])
diss_lo_leakage = qb.get_power(sa, freq=qb.pars['diss_LO'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
diss_im_leakage = qb.get_power(sa, freq=qb.pars['diss_LO']-qb.pars['diss_IF'],span = 1e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
diss_on_power = qb.get_power(sa, freq=qb.pars['diss_LO']+qb.pars['diss_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
qb.opt_mixer( sa, cal='LO',  freq_span = 1e6, reference = ref_L, mode='fine',element='rr')
qb.opt_mixer( sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element='rr')

#%% time of flight
set_attenuator(0)
qb.tof_cal()

#%% resonator spectroscopy

# I,Q,freqs,job = qb.run_scan(df = 25e3, n_avg = 250, element='resonator', chunksize = 90e6, attenuation=35, lo_min = 5.636e9, lo_max = 5.636e9,
#            showprogress=True, res_ringdown_time = int(4e3))
inst.turn_off_ffl_drive()
inst.set_rr_LO(qb.pars['rr_LO'])
qb.resonator_spec(f_LO=5.636e9,atten=35,IF_min=30e6,IF_max=60e6,df=0.1e6,n_avg=1000,savedata=True)

#%% qubit spectroscopy

I,Q,freqs,job = qb.run_scan(df = 50e3, n_avg = 1000, element='qubit', chunksize = 350e6,  lo_min = 3.5e9, lo_max = 3.5e9,
            amp_q_scaling = 0.7, saturation_dur = 20e3, showprogress=True, res_ringdown_time = int(4e3))

I,Q,freqs,job = qb.qubit_spec(f_LO=4.4e9,amp_q_scaling=0.9,IF_min=50e3,IF_max=400e6,df=50e3,n_avg=1000,on_off=False,showprogress=True,savedata=False,check_mixers=False,)
curve_fit_gui(None, freqs, np.abs(I +  1j*Q))

#%% Rabi

qb.update_value('pi_amp', 0.4)
qb.update_value('pi_len', 64)
amps, I, Q, job, period = qb.power_rabi(a_min = 0.01, a_max = 1.0, da = 0.005,  n_avg = 5000,
                                        fit = True, plot = True, detuning = 0e6, check_mixers = False)


t_arr, I, Q, job, pi2width = qb.pulse_exp(exp='rabi',n_avg = 10000, tmin = 0,tmax = 400, dt = 4, amp_q_scaling = 1, fit = True, plot = True,
                                     detuning = 0e6, check_mixers=False)
qb.write_pars()
#%% optimize readout


#%% Ramsey
inst.turn_off_ffl_drive()
inst.turn_on_ffl_drive()
t_arr, I, Q, job, fitted_pars = qb.pulse_exp(exp='ramsey', check_mixers = False, n_avg=5000,dt=4,tmin=0,tmax=2e3,detuning=0e3)

#%% Echo
qb.pulse_exp(exp = 'echo', n_avg = 1000, tmin = qb.pars['pi_half_len'], tmax = 16 * 100+qb.pars['pi_half_len'], dt = 16, fit=True, check_mixers=False)

#%% T1

t_arr, I, Q,job,fitted_pars = qb.pulse_exp(exp='T1',n_avg=5000,tmin=0,tmax=10e3,dt=40,detuning=0)
#%% resonator ring down
qb.pulse_exp(exp = 'ringdown', n_avg = 1000, tmin = 16, tmax =2000, dt = 16, fit=False, check_mixers=False)


#%% diss T1
inst.turn_off_ffl_drive()

qb.update_value('rr_atten', 20)
qb.update_value('rr_pulse_len_in_clk_diss', 100)
t_arr, I, Q, job, fitted_pars = qb.pulse_exp(exp = 'dissT1', n_avg = 20e6, tmin = 16, tmax =200, dt = 4, fit=True, check_mixers=False)

#%% diss spec
I,Q,freqs,job = qb.qubit_spec(f_LO=4.4e9,amp_q_scaling=0.9,IF_min=50e3,IF_max=400e6,df=50e3,n_avg=1000,on_off=False,showprogress=True,savedata=False,check_mixers=False,)

#%% cavity ring down
inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=30e6,IF_max=60e6,df=0.1e6,n_avg=1000,savedata=True)
fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
qb.update_value('rr_freq', fc)
dataDict, fig = measure_ringdown_drive_off(tmax=4e3, n_avg=5000)
dataDict, fig = measure_ringdown_drive_on(amp_ffl_scale=0.05,tmax=1e3, dt=4, n_avg=10000)

#%% qubit reset
dataDict, fig = measure_t1_ffl_off(amp_r_scale=1,amp_ffl_scale=0.3, tmin = 0, tmax = 15e3, dt = 64, n_avg = 5000,)
# qb.update_value('ffl_freq', 3.8918e9)
qb.update_value('ffl_freq', 3.9918e9)
dataDict, fig = measure_t1_w_ffl(amp_r_scale=1,amp_ffl_scale=0.2, tmin = 0, tmax = 8e3, dt = 42, n_avg = 10000,)

# state prep
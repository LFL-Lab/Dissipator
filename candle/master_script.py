# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:42:53 2022

@author: lfl
"""
from mixer_calibration import get_power,opt_mixer
import measurement_functions as mf

with open('pars.json', 'r') as openfile:
    pars = json.load(openfile)
#%% Mixer Calibration
'''Qubit mixer calibration'''
ref_H = - 10
ref_L = -45
mixer1 = 'qubit'
# get LO leakage power
qb_lo_leakage = get_power(sa,freq=qb_LO,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
# get image leakage power
qb_im_leakage = get_power(sa, freq=qb_LO-qb_IF,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
# get qubit drive power at (almost) top of fridge
qb_on_power = get_power(sa, freq=qb_LO+qb_IF,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
opt_mixer(sa, cal='LO',mode='coarse', freq_span = 1e6, reference=ref_H,element=mixer1)
# do a finer sweep in desired range
opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine',element=mixer1)
# do a coarse sweep to minimize sideband
opt_mixer(sa, cal='SB', freq_span = 1e6, reference=ref_H, mode='coarse',element=mixer1)
# do a finer sweep in desired range
opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element=mixer1)

'''Do IQ imbalance calibration'''

mixer2 = 'rr'
# get LO leakage power
rr_lo_leakage = get_power(sa, freq=rr_LO,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
# get image leakage power
rr_im_leakage = get_power(sa, freq=rr_LO-rr_IF,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
# get qubit drive power at (almost) top of fridge
rr_on_power = get_power(sa, freq=rr_LO+rr_IF,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
# do a coarse sweep to minimize LO leakage
opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element=mixer2)
# do a finer sweep in desired range
opt_mixer( sa, cal='LO',  freq_span = 1e6, reference = ref_L, mode='fine',element=mixer2)
# do a coarse sweep to minimize sideband
opt_mixer( sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element=mixer2)
# do a finer sweep in desired range
opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element=mixer2)

sa_close_device(sa)

#%% spectroscopy

#%% Rabi

#%% Ramsey

#%% Echo

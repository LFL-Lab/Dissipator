# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:03:43 2022

Description: A script that generates a JSON file containing mixer calibration parameters
and qubit/resonator frequencies. This JSON file is used by the configuration file.

@author: Evangelos Vlachos
"""

import json
from numpy import pi

# Data to be writen
# qubit mixer calibration parameters
par_dict = {
    "qb_LO":                int(3.8e9),
    "qb_freq":              int(3.32402e9),
    "qb_offsets":           [-0.015517,-0.005172], # I,Q
    "qb_imbalance":         [0.077586,-0.108620], # gain,phase
    "pi_half_len":          260, # needs to be multiple of 4
    "pi_amp":               0.4383,
    "amp_q":                0.45,
    "gauss_len":            48,
    "gauss_amp":            0.45,
    "rr_LO":                int(7.3e9),
    "rr_freq":              int(7.2578e9),
    "rr_offsets":           [-0.0086,-0.00517],
    "rr_imbalance":         [0.015517,0.025862],
    "amp_r":                0.45,
    "rr_pulse_len_in_clk":  2500, # length of readout integration weights in clock cycles 
    "IQ_rotation":          -0/180*pi, # phase rotation applied to IQ data
    "analog_input_offsets": [-0.089918,-0.090406] 
    }

# write dictionary in JSON format without actually creating the file
with open("pars.json", "w") as outfile:
    json.dump(par_dict, outfile)
    
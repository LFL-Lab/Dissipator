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
try:
    print('Opening parameter JSON file')
    with open('pars.json', 'r') as openfile:
        pars = json.load(openfile)
except FileNotFoundError:
        print('JSON file not found; Creating new file using template')
        pars = {
            "qb_LO":                        int(3.0e9),
            "qb_freq":                      int(4e9),
            "qubit_mixer_offsets":          [-0.015517,-0.005172], # I,Q
            "qubit_mixer_imbalance":        [0.077586,-0.108620], # gain,phase
            "pi_half_len":                  260, # needs to be multiple of 4
            "pi_amp":                       0.4383,
            "amp_q":                        0.45,
            "gauss_len":                    48,
            "gauss_amp":                    0.45,
            "rr_LO":                        int(6.55e9),
            "rr_freq":                      int(6.65775e9),
            "rr_mixer_offsets":             [-0.0086,-0.00517],
            "rr_mixer_imbalance":           [0.015517,0.025862],
            "amp_r":                        0.45,
            "rr_pulse_len_in_clk":          2500, # length of readout integration weights in clock cycles
            "IQ_rotation":                  -0/180*pi, # phase rotation applied to IQ data
            "analog_input_offsets":         [-0.089918,-0.090406]
            }

        # write dictionary in JSON format without actually creating the file
        with open("pars.json", "w") as outfile:
            json.dump(pars, outfile)

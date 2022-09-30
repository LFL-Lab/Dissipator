# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:59:17 2022

Description: A script to initialize and setup instruments other than the OPX
@author: Evangelos Vlachos (evlachos@usc.edu)
"""

import numpy as np
from config import *
import sys
sys.path.append("D:\Program Files\Keysight\Labber\Script")
sys.path.append(r"C:\Users\lfl\measurements_test")
import Labber


client = Labber.connectToServer()

# initialize qubit LO 
qLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Qubit', startup = 'Get config'))
qLO.startInstrument()
qLO.setValue('Frequency', qb_LO)
qLO.setValue('Output',True)
config['elements']['qubit']['mixInputs']['lo_frequency'] = qb_LO
config['mixers']['mixer_q1'][0]['lo_frequency'] = qb_LO

# initialize readout LO 
rrLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Readout', startup = 'Get config'))
rrLO.startInstrument()
rrLO.setValue('Frequency', rr_LO)
rrLO.setValue('Output',True)
config['elements']['rr']['mixInputs']['lo_frequency'] = rr_LO
config['mixers']['mixer_rl1'][0]['lo_frequency'] = rr_LO

# initialize digital attenuator 
attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='Readout',address='26777'))
attn.startInstrument()
attn.setValue('Attenuation',17)

# initialize Keithley
# SC = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='24'))
# SC.startInstrument()


if name == '__main__':
    pass

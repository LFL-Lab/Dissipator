# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:59:17 2022

Description: A script to initialize and setup instruments other than the OPX
@author: Evangelos Vlachos (evlachos@usc.edu)
"""

import sys
from VISAdrivers.sa_api import *
sys.path.append("D:\Program Files\Keysight\Labber\Script")
sys.path.append(r"C:\Users\lfl\measurements_test")
import Labber
import json

client = Labber.connectToServer()
# qmm = QuantumMachinesManager()
# qm = qmm.open_qm(config)

with open('pars.json', 'r') as openfile:
    pars = json.load(openfile)

def set_qb_LO(freq):
    qb_LO = pars['qb_LO']
    print(f'Setting qubit LO to {round(freq*1e-9,5)} GHz')
    # initialize qubit LO
    qLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Qubit', startup = 'Get config'))
    qLO.startInstrument()
    qLO.setValue('Frequency', freq)
    qLO.setValue('Output',True)

def set_rr_LO(freq):
    # rr_LO = pars['rr_LO']
    print(f'Setting readout LO to {round(freq*1e-9,5)} GHz')
    # initialize readout LO
    rrLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Readout', startup = 'Get config'))
    rrLO.startInstrument()
    rrLO.setValue('Frequency', freq)
    rrLO.setValue('Output',True)

def set_attenuator(attenuation):
    # initialize digital attenuator
    attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='Readout',address='26777'))
    attn.startInstrument()
    attn.setValue('Attenuation',attenuation)

def get_attenuation():
    # initialize digital attenuator
    attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='Readout',address='26777'))
    attn.startInstrument()

    return attn.getValue('Attenuation')

def init_sa():
    sa = sa_open_device()["handle"]
    sa_config_acquisition(device = sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
    sa_config_gain_atten(sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
    sa_config_sweep_coupling(device = sa, rbw = 1e3, vbw = 1e3, reject=0)

    return sa
# initialize Keithley
# SC = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='24'))
# SC.startInstrument()

if __name__ == '__main__':
    # initialize spectrum analyzer
    set_qb_LO(freq=pars['qb_LO'])
    set_rr_LO(freq=pars['rr_LO'])
    set_attenuator(0)
    sa = init_sa()
else:
    pass

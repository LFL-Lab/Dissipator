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

client = Labber.connectToServer()
qmm = QuantumMachinesManager()
qm = qmm.open_qm(config)

with open('pars.json', 'r') as openfile:
    pars = json.load(openfile)

qb_LO = pars['qb_LO']
# initialize qubit LO
qLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Qubit', startup = 'Get config'))
qLO.startInstrument()
qLO.setValue('Frequency', qb_LO)
qLO.setValue('Output',True)


rr_LO = pars['rr_LO']
# initialize readout LO
rrLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Readout', startup = 'Get config'))
rrLO.startInstrument()
rrLO.setValue('Frequency', rr_LO)
rrLO.setValue('Output',True)

# initialize digital attenuator
attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='Readout',address='26777'))
attn.startInstrument()
attn.setValue('Attenuation',17)

# initialize spectrum analyzer
sa = sa_open_device()["handle"]
sa_config_acquisition(device = sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
sa_config_gain_atten(sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
sa_config_sweep_coupling(device = sa, rbw = 1e3, vbw = 1e3, reject=0)

# initialize Keithley
# SC = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='24'))
# SC.startInstrument()


if __name__ == '__main__':
    pass

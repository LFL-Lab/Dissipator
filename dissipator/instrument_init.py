# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:59:17 2022

Description: A script to initialize and setup instruments other than the OPX
@author: Evangelos Vlachos (evlachos@usc.edu)
"""
# import os
# os.chdir()
import sys
from VISAdrivers.sa_api import *
sys.path.append("D:\Program Files\Keysight\Labber\Script")
sys.path.append(r"C:\Users\lfl\measurements_test")
import Labber
import json

client = Labber.connectToServer()
# qmm = QuantumMachinesManager()
# qm = qmm.open_qm(config)
qubit_LO_model = 'SignalCore SC5506A Signal Generator'
qubit_LO_name = '10002A08' 
qubit_LO_quantity_name = {'freq':'RF1 frequency',
                          'power': 'RF1 power level',
                          'output': 'RF1 output status'}

rr_LO_model = 'SignalCore SC5511A Signal Generator'
rr_LO_name = '10002A05'
rr_LO_quantity_name = {'freq': 'Frequency',
                       'power': 'Power',
                       'output':'Output status'}

ffl_LO_model = 'BNC 845 Signal Generator'
ffl_LO_name = 'FFL drive'
ffl_LO_quantity_name = {'freq': 'Frequency',
                       'output':'Output'}

diss_LO_model = 'SignalCore SC5511A Signal Generator'
diss_LO_name = '10002F25'
diss_LO_quantity_name = {'freq': 'Frequency',
                        'power': 'Power',
                        'output':'Output status'}

def set_qb_LO(freq):
    print(f'Setting qubit LO to {round(freq*1e-9,5)} GHz')
    # initialize qubit LO
    qLO = client.connectToInstrument(qubit_LO_model, dict(name=qubit_LO_name, startup = 'Get config'))
    qLO.startInstrument()
    qLO.setValue(qubit_LO_quantity_name['freq'], freq)
    qLO.setValue(qubit_LO_quantity_name['power'], 15)
    qLO.setValue(qubit_LO_quantity_name['output'],True)

def get_qb_LO():
    # initialize qubit LO
    qLO = client.connectToInstrument(qubit_LO_model, dict(name=qubit_LO_name, startup = 'Get config'))
    qLO.startInstrument()
    return qLO.getValue(qubit_LO_quantity_name['freq'])

def set_diss_LO(freq):
    print(f'Setting dissipator LO to {round(freq*1e-9,5)} GHz')
    # initialize qubit LO
    dLO = client.connectToInstrument(diss_LO_model, dict(name=diss_LO_name, startup = 'Get config'))
    dLO.startInstrument()
    dLO.setValue(diss_LO_quantity_name['freq'], freq)
    dLO.setValue(diss_LO_quantity_name['power'], 15)
    dLO.setValue(diss_LO_quantity_name['output'],True)

def get_diss_LO():
    # initialize qubit LO
    dLO = client.connectToInstrument(diss_LO_model, dict(name=diss_LO_name, startup = 'Get config'))
    dLO.startInstrument()
    return qLO.getValue(diss_LO_quantity_name['freq'])

def set_rr_LO(freq):
    print(f'Setting readout LO to {round(freq*1e-9,5)} GHz')
    # initialize readout LO
    rrLO = client.connectToInstrument(rr_LO_model, dict(name=rr_LO_name, startup = 'Get config'))
    rrLO.startInstrument()
    rrLO.setValue('Frequency', freq)
    rrLO.setValue(rr_LO_quantity_name['power'], 15)
    rrLO.setValue(rr_LO_quantity_name['output'],True)

def get_rr_LO():
    # initialize qubit LO
    rrLO = client.connectToInstrument(rr_LO_model, dict(name=rr_LO_name, startup = 'Get config'))
    rrLO.startInstrument()
    return rrLO.getValue('Frequency')

def get_ffl_LO():
    # initialize qubit LO
    fLO = client.connectToInstrument(ffl_LO_model, dict(name=ffl_LO_name, startup = 'Get config'))
    fLO.startInstrument()
    return fLO.getValue(ffl_LO_quantity_name['freq'])

def set_ffl_LO(freq):
    print(f'Setting ffl LO to {round(freq*1e-9,5)} GHz')
    # initialize ffl LO
    fLO = client.connectToInstrument(ffl_LO_model, dict(name=ffl_LO_name, startup = 'Get config'))
    fLO.startInstrument()
    fLO.setValue('Frequency', freq)
    fLO.setValue(ffl_LO_quantity_name['output'],True)
    

    
def turn_on_ffl_drive():
    fLO = client.connectToInstrument(ffl_LO_model, dict(name=ffl_LO_name, startup = 'Get config'))
    fLO.startInstrument()
    fLO.setValue(ffl_LO_quantity_name['output'],True)
    
def turn_off_ffl_drive():
    fLO = client.connectToInstrument(ffl_LO_model, dict(name=ffl_LO_name, startup = 'Get config'))
    fLO.startInstrument()
    fLO.setValue(ffl_LO_quantity_name['output'],False)
    
def turn_diss_LO_off():
    dLO = client.connectToInstrument(diss_LO_model, dict(name=diss_LO_name, startup = 'Get config'))
    dLO.startInstrument()
    dLO.setValue(diss_LO_quantity_name['output'],False)
    
def turn_qb_LO_off():
    qLO = client.connectToInstrument(qubit_LO_model, dict(name=qubit_LO_name, startup = 'Get config'))
    qLO.startInstrument()
    qLO.setValue(qubit_LO_quantity_name['output'],False)
    
def set_attenuator(attenuation):
    # initialize digital attenuator
    attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='readout attenuator',address='26776'))
    attn.startInstrument()
    attn.setValue('Attenuation',attenuation)

def get_attenuation():
    # initialize digital attenuator
    attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='readout attenuator',address='26776'))
    attn.startInstrument()

    return attn.getValue('Attenuation')

def init_sa():
    sa = sa_open_device()["handle"]
    sa_config_acquisition(device = sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
    sa_config_gain_atten(sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
    sa_config_sweep_coupling(device = sa, rbw = 1e2, vbw = 1e2, reject=0)
    return sa

def init_sa_by_serial_number(serial_number=20234492):
    sa = sa_open_device_by_serial(serial_number)["handle"]
    sa_config_acquisition(device = sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
    sa_config_gain_atten(sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
    sa_config_sweep_coupling(device = sa, rbw = 1e2, vbw = 1e2, reject=0)
    return sa

# initialize Keithley
# SC = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='24'))
# SC.startInstrument()

# try:
#     sa
# except NameError:
#     sa = init_sa()

if __name__ == '__main__':
    # initialize spectrum analyzer
    sa = init_sa()
    # sa2 = init_sa_by_serial_number()
else:
    pass
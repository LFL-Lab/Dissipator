# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 10:29:59 2021

@author: lfl
"""
from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from config import *
import matplotlib.pyplot as plt
import numpy as np


import sys
sys.path.append("D:\Program Files\Keysight\Labber\Script")
import Labber




client = Labber.connectToServer()

qLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Qubit', startup = 'Get config'))
qLO.startInstrument()
qLO.setValue('Frequency', qb_LO)
qLO.setValue('Output',True)
config['elements']['qubit']['mixInputs']['lo_frequency'] = qb_LO
config['mixers']['mixer_q1'][0]['lo_frequency'] = qb_LO

# rrLO = client.connectToInstrument('SignalCore SC5511A Signal Generator', dict(name='10002F1D', startup = 'Get config'))
rrLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='Readout', startup = 'Get config'))
rrLO.startInstrument()
rrLO.setValue('Frequency', rr_LO)
rrLO.setValue('Output',True)
config['elements']['rr']['mixInputs']['lo_frequency'] = rr_LO
config['mixers']['mixer_rl1'][0]['lo_frequency'] = rr_LO


###################
# The QUA program #
###################


with program() as play_pulses:
    with infinite_loop_():
        play("const", 'qubit',duration=100)
        play("readout", "rr", duration=100)

    

######################################
# Open Communication with the Server #
######################################
qmm = QuantumMachinesManager()

####################
# Simulate Program #
####################
# simulation_config = SimulationConfig(
#                     duration=90000,
#                     simulation_interface=LoopbackInterface([("con1", 3, "con1", 1), ("con1", 4, "con1", 2)]))

# job = qmm.simulate(config, rr_spec, simulation_config)
qm = qmm.open_qm(config)
job = qm.execute(play_pulses)



### Some functions for setting mixer corrections a little more easily

def getmixer(element: str): return config["elements"][element]["mixInputs"]["mixer"]

""" RF input : "I" or "Q" """
def get_mixer_input(element: str, RF_input: str): return config["elements"][element]["mixInputs"][RF_input]

def get_dc_offset(element: str, RF_input: str): 
    
    controller, port = get_mixer_input(element, RF_input)
    return config["controllers"][controller]["analog_outputs"][port]["offset"]

def set_dc_offset(element: str, RF_input: str, value: float):
    qm.set_output_dc_offset_by_element(element, RF_input, value)
    controller, port = get_mixer_input(element, RF_input)
    config["controllers"][controller]["analog_outputs"][port]["offset"] = value


def step_dc_offset(element: str, RF_input: str, step: float):
    offset = get_dc_offset(element, RF_input)
    set_dc_offset(element, RF_input, offset + step)
    
        
def set_IQ_imbalance(element: str, index: int, value: float):
    
    if element == "qubit":
        qubit_imbalance[index] = value
        qm.set_mixer_correction(getmixer(element), qb_IF, qb_LO, IQ_imbalance(*qubit_imbalance))
    elif element == "rr":
        rr_imbalance[index] = value
        qm.set_mixer_correction(getmixer(element), rr_IF, rr_LO, IQ_imbalance(*rr_imbalance))
    else:
        raise Exception("Please enter 'qubit' or 'rr' as the element.")
        

def step_IQ_imbalance(element: str, index: int, step: float):
    
    if element == "qubit":
        qubit_imbalance[index] += step
        qm.set_mixer_correction(getmixer(element), qb_IF, qb_LO, IQ_imbalance(*qubit_imbalance))
    
    elif element == "rr":
        rr_imbalance[index] += step
        qm.set_mixer_correction(getmixer(element), rr_IF, rr_LO, IQ_imbalance(*rr_imbalance))
        
    else:
        raise Exception("Please enter 'qubit' or 'rr' as the element.")
    
            

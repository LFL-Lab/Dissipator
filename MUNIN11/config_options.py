
########
# CONSTANTS 
########

calibration_dict = {
    
    'R1' : {    'res'   :   {   'freq'      : int(6.816860e9),
                                'LO'        : int(6.75e9)
                                'imbalance' : [-0.0156, -0.055]
                                'offset'    : {'I' : -0.0225, 'Q' : 0.006}
                            }
            
                'qubits' :   {  '3.44787e9' : 
                                                  { 'freq'      : int(3.44787e9),
                                                    'LO'        : int(3.40e9),
                                                    'imbalance' : [0.075, -0.12],
                                                    'offset'    : {'I' : -0.0178, 'Q' : -0.0216}
                                                  }
                    
                             }
                    
                            
            
        
        }
    
    
    
    
    
    
    }



res1 = int(6.816860e9)
res2 = int(7.026440e9)
res3 = int(7.105125e9)
res3_biased = int(7.1069e9)

qb_A = int(5.11275e9) 
qb_B = int(5.3176e9)
qubit1 = int(3.44787e9)
qubit2 = int(4.49e9) #int(4.519e9)
qubit3 = int(4.995e9)
qubit3_biased = int(6.233e9) #int(6.2330e9)
qubit3_2 = int(4.78e9)
qubit3_1 = int(4.55e9)
qubit554 = int(5.54e9)

testfreq = int(5.45e9)



########
# SETTINGS
########


#####
# resonator settings
#####
resFreq = res2

if resFreq is res1:
    
    rr_LO = int(6.75e9) # don't forget to change this on the signal generator
    rr_imbalance = [-0.0156, -0.055]
    rrI_analog_out = {"offset": -0.0225} # DC offset on I #-0.025
    rrQ_analog_out = {"offset": 0.006} # DC offset on Q
    res_mixer_IQ_imbalance = IQ_imbalance(*rr_imbalance) 
    
elif resFreq is res2:
    
    rr_LO = int(7.08e9) # don't forget to change this on the signal generator
    rrI_analog_out = {"offset": -0.01} # DC offset on I
    rrQ_analog_out = {"offset": -0.001} # DC offset on Q
    res_mixer_IQ_imbalance = IQ_imbalance(-0.005,-0.06)
    
elif resFreq is res3:
    
    rr_LO = int(7.16e9) # don't forget to change this on the signal generator
    rrI_analog_out = {"offset": -0.02} # DC offset on I
    rrQ_analog_out = {"offset": -0.0035} # DC offset on Q
    res_mixer_IQ_imbalance = IQ_imbalance(-0.005,-0.06)
    
elif resFreq is res3_biased: 
    
    rr_LO = int(7.16e9) # don't forget to change this on the signal generator
    rrI_analog_out = {"offset": -0.02} # DC offset on I
    rrQ_analog_out = {"offset": -0.0035} # DC offset on Q
    res_mixer_IQ_imbalance = IQ_imbalance(-0.005,-0.06)
    

rr_IF = resFreq - rr_LO
amp_r = 0.45 # amplitude of const and readout wfms in V, max is 0.5, but we shouldn't use the full range # was 0.13


#####
# qubit settings
#####
qbFreq = qubit2

if qbFreq is qb_A:
    
    qb_LO = int(5.05e9) # don't forget to change this on the signal generator
    qb_mixer_IQ_imbalance = IQ_imbalance(-0.02,-0.05) # need to calibrate; remove note when calibrated
    qbI_analog_out = {"offset": -0.018} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": 0.0001} # need to calibrate; remove note when calibrated

elif qbFreq is qb_B:
    
    qb_LO = int(5.287e9) # don't forget to change this on the signal generator
    qb_mixer_IQ_imbalance = IQ_imbalance(0.052,0.006) # need to calibrate
    qbI_analog_out = {"offset": -0.01} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": -0.0004} # need to calibrate; remove note when calibrated 
    
elif qbFreq is qubit1:
    
    qb_LO = int(3.40e9)  # don't forget to change this on the signal generator
    qubit_imbalance = [0.075, -0.12]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.0178} 
    qbQ_analog_out = {"offset": -0.0216} 
    
    pi_half_len = 76 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =   0.41 #0.274
    
    
elif qbFreq is qubit2:
    
    # lower sideband
    # qb_LO = int(4.60e9)
    # qubit_imbalance = [0.041, 0.098]
    # qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    # qbI_analog_out = {"offset": -0.01} # need to calibrate; remove note when calibrated
    # qbQ_analog_out = {"offset": -0.0084} # need to calibrate; remove note when calibrated
    
    
    # upper sideband
    qb_LO = int(4.45e9)  # don't forget to change this on the signal generator
    qubit_imbalance = [0.096, 0.0]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.0189} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": -0.0061} # need to calibrate; remove note when calibrated
    
    pi_half_len = 32 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =   0.3597 #0.274
    
    
elif qbFreq is qubit3_biased:
    
    # upper sideband
    qb_LO = int(6.0e9)  # don't forget to change this on the signal generator
    qubit_imbalance = [0.026, -0.022]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.016} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": -0.0042} # need to calibrate; remove note when calibrated
        
    pi_half_len = 32 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =   0.3597 #0.274
    
    
elif qbFreq is qubit3:
    
    # for lower sideband generation
    qb_LO = int(5.05e9)
    qb_mixer_IQ_imbalance = IQ_imbalance(0.0718, 0.019)
    qbI_analog_out = {"offset": -0.0081} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": -0.0014} # need to calibrate; remove note when calibrated
    
    # for upper sideband generation
    # qb_LO = int(4.92e9)  # don't forget to change this on the signal generator
    # qb_mixer_IQ_imbalance = IQ_imbalance(0.067,0.038)
    # qbI_analog_out = {"offset": -0.007} # need to calibrate; remove note when calibrated
    # qbQ_analog_out = {"offset": -0.0021} # need to calibrate; remove note when calibrated
    
elif qbFreq is qubit3_2:
    
    qb_LO = int(4.73e9)  # don't forget to change this on the signal generator
    qb_mixer_IQ_imbalance = IQ_imbalance(0.067,0.038)
    qbI_analog_out = {"offset": -0.007} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": -0.0021} # need to calibrate; remove note when calibrated

elif qbFreq is qubit3_1:
    
    qb_LO = int(4.50e9)  # don't forget to change this on the signal generator
    qb_mixer_IQ_imbalance = IQ_imbalance(0.071,0.067)
    qbI_analog_out = {"offset": -0.008} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": -0.0076} # need to calibrate; remove note when calibrated
    
elif qbFreq is qubit554:
    
    qb_LO = int(5.48e9)  # don't forget to change this on the signal generator
    qb_mixer_IQ_imbalance = IQ_imbalance(0.075,-0.015)
    qbI_analog_out = {"offset": -0.011} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": 0.0011} # need to calibrate; remove note when calibrated
    
    pi_half_len = 76 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =   0.41 #0.274
  
    
    
elif qbFreq is testfreq:
    
    qb_LO = int(5.50e9)  # don't forget to change this on the signal generator
    qubit_imbalance = [0.0635, 0.0215]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.0108} # need to calibrate; remove note when calibrated
    qbQ_analog_out = {"offset": 0.0} # need to calibrate; remove note when calibrated
    

qb_IF = qbFreq - qb_LO
amp_q = 0.45 # 0.45 is the maximum, can make smaller to see if resonances change shape





#####
# other settings
#####

rr_pulse_len_in_clk = 500 # needs to be a multiple of 4

phi = -0/180 * np.pi # the angle to rotate the IQ plane such that most of the information lies on a single quadrature

# pi pulse depends on the qubit, so moved into qubit-dependent code above
# eventually would be good to handle this using a dictionary

# pi_half_len = 32 # needs to be a multiple of 4, in ns unit
# pi_len = 2 * pi_half_len
# pi_amp =   0.3597 #0.274

gauss_len = 48
gauss_amp = 0.44 #0.0665
gauss_half_amp = gauss_amp / 2

gauss_wf_4ns = delayed_gauss(gauss_amp, 4, 2)





#####
# configuration dict for QUA interpretor
#####
config = {

    "version": 1,

    "controllers": {
        "con1": {
            "type": "opx1",
            "analog_outputs": {
                1: rrI_analog_out,  # rr I
                2: rrQ_analog_out,  # rr Q  
                3: qbI_analog_out,  # qubit I
                4: qbQ_analog_out,  # qubit Q
            },
            "digital_outputs": {},
            "analog_inputs": {
                1: {"offset": -0.09739892578125 - 0.002783935546875, "gain_db": 3},  # rr I
                2: {"offset": -0.0818253173828125 + 0.0008167724609375, "gain_db": 3}  # rr Q
            },
        },
    },
    
    "elements": {
        "qubit": {
            "mixInputs": {
                "I": ("con1", 3),
                "Q": ("con1", 4),
                "lo_frequency": qb_LO,
                "mixer": "mixer_q1",
            },
            "intermediate_frequency": qb_IF,
            "digitalInputs": {},
            "operations": {
                "const": "const_pulse_IQ",
                "gauss": "gaussian_pulse",
                "gauss_4ns": "gaussian_4ns_pulse",
                "pi": "pi_pulse1",
                "pi_half": "pi_half_pulse1",
                "arb_op": "arb_pulse",
                "X": "Xpi_pulse",
                "Y": "Ypi_pulse",
                "X/2": "Xpi_half_pulse",
                "Y/2": "Ypi_half_pulse",
            },
        },
        "rr": {
            "mixInputs": {
                "I": ("con1", 1),
                "Q": ("con1", 2),
                "lo_frequency": rr_LO,
                "mixer": "mixer_rl1",
            },
            "intermediate_frequency": rr_IF,
            "outputs": {
                "out1": ("con1", 1),
                "out2": ("con1", 2),
            },
            "time_of_flight": 276, #288, # should be multiple of 4 (at least 24)
            "smearing": 0, # adds 40ns of data from each side of raw adc trace to account for ramp up and down of readout pulse
            "operations": {
                "const": "const_pulse_IQ_rr",
                "readout": "ro_pulse1",
            },
        },
    },
    
    "pulses": {
        "const_pulse_IQ": {
            "operation": "control",
            "length": 100,
            "waveforms": {
                "I": "const_wf",
                "Q": "zero_wf",
            },
        },
        "const_pulse_IQ_rr": {
            "operation": "control",
            "length": 100,
            "waveforms": {
                "I": "const_wf_rr",
                "Q": "zero_wf",
            },
        },
        "pi_pulse1": {
            "operation": "control",
            "length": pi_len,
            "waveforms": {
                "I": "pi_wf_i1",
                "Q": "pi_wf_q1",
            },
        },
        "Xpi_pulse": {
            "operation": "control",
            "length": pi_len,
            "waveforms": {
                "I": "pi_wf_i1",
                "Q": "pi_wf_q1",
            },
        },
        "Ypi_pulse": {
            "operation": "control",
            "length": pi_len,
            "waveforms": {
                "I": "pi_wf_q1",
                "Q": "pi_wf_i1",
            },
        },
        "Xpi_half_pulse": {
            "operation": "control",
            "length": pi_half_len,
            "waveforms": {
                "I": "pi_half_wf_i1",
                "Q": "pi_half_wf_q1",
            },
        },
        "Ypi_half_pulse": {
            "operation": "control",
            "length": pi_half_len,
            "waveforms": {
                "I": "pi_half_wf_q1",
                "Q": "pi_half_wf_i1",
            },
        },
        "gaussian_pulse": {
            "operation": "control",
            "length": gauss_len,
            "waveforms": {
                "I": "gaussian_wf",
                "Q": "zero_wf",
            },
        },
        "gaussian_4ns_pulse": {
            "operation": "control",
            "length": 16,
            "waveforms": {
                "I": "gaussian_4ns_wf",
                "Q": "zero_wf",
            },
        },
        "pi_half_pulse1": {
            "operation": "control",
            "length": pi_half_len,
            "waveforms": {
                "I": "pi_half_wf_i1",
                "Q": "pi_half_wf_q1",
            },
        },
        "ro_pulse1": {
            "operation": "measurement",
            "length": rr_pulse_len_in_clk *4, # in ns (needs to be multiple of 4)
            "waveforms": {"I": "ro_wf1", "Q": "zero_wf"},
            "integration_weights": {
                "integW_cos": "integW1_cos",
                "integW_sin": "integW1_sin",
                "integW_minus_sin": "integW1_minus_sin"
            },
            "digital_marker": "ON",
        },
        "arb_pulse": {
            "operation": "control",
            "length": 40,
            "waveforms": {
                "I": "arb_wfm",
                "Q": "zero_wf",
            },
        },
    },

    "waveforms": {
        "zero_wf": {"type": "constant", "sample": 0.0},
        "const_wf": {"type": "constant", "sample": amp_q},
        "const_wf_rr": {"type": "constant", "sample": amp_r},
        "gaussian_wf": {"type": "arbitrary", "samples": [float(arg) for arg in gauss_amp * gaussian(gauss_len, gauss_len/5)]},
        "gaussian_4ns_wf": {"type": "arbitrary", "samples": gauss_wf_4ns},
        "ro_wf1": {"type": "constant", "sample": amp_r},
        "pi_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pi_amp * gaussian(pi_len, pi_len/5)]},
        "pi_wf_q1": {"type": "constant", "sample": 0.0},
        "pi_half_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pi_amp * gaussian(pi_half_len, pi_half_len/5)]},
        "pi_half_wf_q1": {"type": "constant", "sample": 0.0},
        "arb_wfm": {"type": "arbitrary", "samples": [0.2]*10+[0.3]*10+[0.25]*20},
    },

    "digital_waveforms": {
        "ON": {"samples": [(1, 0)]}
    },

    "integration_weights": {
        "integW1_cos": {
            "cosine": [(np.cos(phi) , rr_pulse_len_in_clk*4)],
            "sine": [(-np.sin(phi) , rr_pulse_len_in_clk*4)],
        },
        "integW1_sin": {
            "cosine": [(np.sin(phi) , rr_pulse_len_in_clk*4)],
            "sine": [(np.cos(phi) , rr_pulse_len_in_clk*4)],
        },
        "integW1_minus_sin": {
            "cosine": [(-np.sin(phi) , rr_pulse_len_in_clk*4)],
            "sine": [(-np.cos(phi) , rr_pulse_len_in_clk*4)],
        },
        "integW2_cos": {
            "cosine": [1.0] * rr_pulse_len_in_clk,
            "sine": [0.0] * rr_pulse_len_in_clk,
        },
        "integW2_sin": {
            "cosine": [0.0] * rr_pulse_len_in_clk,
            "sine": [1.0] * rr_pulse_len_in_clk,
        },
        "integW2_minus_sin": {
            "cosine": [0.0] * rr_pulse_len_in_clk,
            "sine": [-1.0] * rr_pulse_len_in_clk,
        }
    },

    "mixers": {
        "mixer_q1": [{"intermediate_frequency": qb_IF, "lo_frequency": qb_LO, "correction": qb_mixer_IQ_imbalance}],
        "mixer_rl1": [{"intermediate_frequency": rr_IF, "lo_frequency": rr_LO, "correction": res_mixer_IQ_imbalance}],
    }
}

# qmm = QuantumMachinesManager()
# qm = qmm.open_qm(config)
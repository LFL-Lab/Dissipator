from scipy.signal.windows import gaussian
import numpy as np
from datetime import date
# from qm.QuantumMachinesManager import QuantumMachinesManager

today = date.today()
datestr = today.strftime("%Y%m%d")

def IQ_imbalance(g, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1-g**2)*(2*c**2-1))
    return [float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]

def delayed_gauss(amp, length, sigma):
    gauss_arg = np.linspace(-sigma, sigma, length)
    delay = 16 - length - 4
    if delay < 0:
        return amp * np.exp(-(gauss_arg ** 2) / 2)

    return np.r_[np.zeros(delay), amp * np.exp(-(gauss_arg ** 2) / 2), np.zeros(4)]

########
# CONSTANTS 
########

R5 = int(7.2578e9)
Q5 = int(3.32402e9)

R4 = int(7.22496e9)
Q4 = int(3.2075e9)

R3 = int(7.2019e9)
Q3 = int(3.09244e9)

R2 = int(7.1726e9)
Q2 = int(3.03361e9)
Q22 = int(3.056e9)

R1 = int(7.13135e9)
Q1 = int(2.81718e9)

R0 = int(7.11885e9)
Q0 = int(2.8034e9)

resFreq = R5
qbFreq = Q5



########
# SETTINGS
########


#####
# resonator settings
#####
rr_LO = int(7.30e9) #int(7.2578e9) # don't forget to change this on the signal generator
rr_imbalance = [0.015517241379310348, 0.025862068965517238]
rrI_analog_out = {"offset": -0.0086} # DC offset on I #-0.025
rrQ_analog_out = {"offset":  -0.00517} # DC offset on Q
res_mixer_IQ_imbalance = IQ_imbalance(*rr_imbalance)

rr_IF = resFreq - rr_LO
amp_r = 0.45 # amplitude of const and readout wfms in V, max is 0.5, but we shouldn't use the full range # was 0.13


#####
# qubit settings
#####
if  qbFreq is Q5:
    qb_LO = int(3.38e9)  # don't forget to change this on the signal generator

    # upper sideband
    qubit_imbalance = [0.07758620689655174, -0.10862068965517241]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.015517241379310341} 
    qbQ_analog_out = {"offset": -0.005172413793103445} 
    
    # lower sideband
    # qubit_imbalance = [0.0569, -0.06724]
    # qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    # qbI_analog_out = {"offset": -0.01897} 
    # qbQ_analog_out = {"offset": -0.01552} 
    
    pi_half_len = 260 #192 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp = 0.4383 #0.274
    
    qb_IF = qbFreq - qb_LO
    amp_q = 0.45 # 0.45 is the maximum, can make smaller to see if resonances change shape
    
if qbFreq is Q4:
    qb_LO = int(3.31e9)  # don't forget to change this on the signal generator

    # upper sideband
    qubit_imbalance = [0.0776, -0.0776]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.01552} 
    qbQ_analog_out = {"offset": -0.01552} 
    
    # lower sideband
    # qubit_imbalance = [0.0569, -0.06724]
    # qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    # qbI_analog_out = {"offset": -0.01897} 
    # qbQ_analog_out = {"offset": -0.01552} 
    
    pi_half_len = 180 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =  0.44 #0.274
    
    qb_IF = qbFreq - qb_LO
    amp_q = 0.45 # 0.45 is the maximum, can make smaller to see if resonances change shape

if qbFreq is Q3:
    qb_LO = int(3.31e9)  # don't forget to change this on the signal generator

    # upper sideband
    qubit_imbalance = [0.0776, -0.0776]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.01552} 
    qbQ_analog_out = {"offset": -0.01552} 
    
    # lower sideband
    # qubit_imbalance = [0.0569, -0.06724]
    # qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    # qbI_analog_out = {"offset": -0.01897} 
    # qbQ_analog_out = {"offset": -0.01552} 
    
    pi_half_len = 236 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =  0.44 #0.274
    
    qb_IF = qbFreq - qb_LO
    amp_q = 0.45 # 0.45 is the maximum, can make smaller to see if resonances change shape
    
    resetTime = 400e3 #JL for IQ blobs
    
if qbFreq is Q22:
    qb_LO = int(3.31e9)  # don't forget to change this on the signal generator

    # upper sideband
    qubit_imbalance = [0.0776, -0.0776]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.01552} 
    qbQ_analog_out = {"offset": -0.01552} 
    
    # lower sideband
    # qubit_imbalance = [0.0569, -0.06724]
    # qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    # qbI_analog_out = {"offset": -0.01897} 
    # qbQ_analog_out = {"offset": -0.01552} 
    
    pi_half_len = 232 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =  0.44 #0.274
    
    qb_IF = qbFreq - qb_LO
    amp_q = 0.45 # 0.45 is the maximum, can make smaller to see if resonances change shape
    
if qbFreq is Q1:
    qb_LO = int(3.31e9)  # don't forget to change this on the signal generator

    # upper sideband
    qubit_imbalance = [0.0776, -0.0776]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.01552} 
    qbQ_analog_out = {"offset": -0.01552} 
    
    # lower sideband
    # qubit_imbalance = [0.0569, -0.06724]
    # qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    # qbI_analog_out = {"offset": -0.01897} 
    # qbQ_analog_out = {"offset": -0.01552} 
    
    pi_half_len = 1392 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =  0.44 #0.274
    
    qb_IF = qbFreq - qb_LO
    amp_q = 0.45 # 0.45 is the maximum, can make smaller to see if resonances change shape    

if qbFreq is Q0:
    qb_LO = int(3.31e9)  # don't forget to change this on the signal generator

    # upper sideband
    qubit_imbalance = [0.0776, -0.0776]
    qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    qbI_analog_out = {"offset": -0.01552} 
    qbQ_analog_out = {"offset": -0.01552} 
    
    # lower sideband
    # qubit_imbalance = [0.0569, -0.06724]
    # qb_mixer_IQ_imbalance = IQ_imbalance(*qubit_imbalance)
    # qbI_analog_out = {"offset": -0.01897} 
    # qbQ_analog_out = {"offset": -0.01552} 
    
    pi_half_len = 1392 # needs to be a multiple of 4, in ns unit
    pi_len = 2 * pi_half_len
    pi_amp =  0.44 #0.274
    
    qb_IF = qbFreq - qb_LO
    amp_q = 0.45 # 0.45 is the maximum, can make smaller to see if resonances change shape  
#####
# other settings
#####

rr_pulse_len_in_clk = 2500 # needs to be a multiple of 4

phi = -0/180 * np.pi # the angle to rotate the IQ plane such that most of the information lies on a single quadrature

# pi pulse depends on the qubit, so moved into qubit-dependent code above
# eventually would be good to handle this using a dictionary

# pi_half_len = 32 # needs to be a multiple of 4, in ns unit
# pi_len = 2 * pi_half_len
# pi_amp =   0.3597 #0.274

gauss_len = 48
# gauss_amp =  0.115 * 0.44
gauss_amp = 0.44
gauss_half_amp = gauss_amp / 2

gauss_wf_4ns = delayed_gauss(gauss_amp, 4, 2)


def res_demod(I, Q):
    
    return (dual_demod.full("integW_cos", "out1", "integW_minus_sin", "out2", I),
            dual_demod.full("integW_sin", "out1", "integW_cos", "out2", Q))


#####
# configuration dict for QUA interpretor
#####
config = {

    "version": 1,

    "controllers": {
        "con1": {
            "type": "opx1",
            "analog_outputs": {
                1: qbI_analog_out,  # qubit I
                2: qbQ_analog_out,  # qubit Q
                3: rrI_analog_out,  # rr I
                4: rrQ_analog_out,  # rr Q  
            },
            "digital_outputs": {},
            "analog_inputs": {
                1: {"offset": -0.08991845070800782, "gain_db": 3},  # rr I
                2: {"offset": -0.09040683835449219, "gain_db": 3}  # rr Q
            },
        },
    },
    
    "elements": {
        "qubit": {
            "mixInputs": {
                "I": ("con1", 1),
                "Q": ("con1", 2),
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
                "I": ("con1", 3),
                "Q": ("con1", 4),
                "lo_frequency": rr_LO,
                "mixer": "mixer_rl1",
            },
            "intermediate_frequency": rr_IF,
            "outputs": {
                "out1": ("con1", 1),
                "out2": ("con1", 2),
            },
            "time_of_flight": 224, # should be multiple of 4 (at least 24)
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
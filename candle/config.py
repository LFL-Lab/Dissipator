from scipy.signal.windows import gaussian
import numpy as np
from datetime import date
# from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
import json
today = date.today()
from utilities import *
datestr = today.strftime("%Y%m%d")

# Opening JSON file containing values of experimental parameters such as optimal DC offsets
with open('pars.json', 'r') as openfile:
    pars = json.load(openfile)

qbFreq = pars['qb_freq']
rrFreq = pars['rr_freq']
qb_LO = pars['qb_LO']
rr_LO = pars['rr_LO']
phi = pars['IQ_rotation']
rr_pulse_len_in_clk = pars['rr_pulse_len_in_clk']
qb_IF = qbFreq - qb_LO
rr_IF = rrFreq - rr_LO
qb_offset_I = {"offset": pars['qubit_mixer_offsets'][0]}
qb_offset_Q = {"offset": pars['qubit_mixer_offsets'][1]}
rr_offset_I = {"offset": pars['rr_mixer_offsets'][0]}
rr_offset_Q = {"offset": pars['rr_mixer_offsets'][1]}

gauss_wf_4ns = delayed_gauss(pars['gauss_amp'], 4, 2)

#####
# configuration dict for QUA interpretor
#####

config = {

    "version": 1,

    "controllers": {
        "con1": {
            "type": "opx1",
            "analog_outputs": {
                1: qb_offset_I,  # qubit I
                2: qb_offset_Q,  # qubit Q
                3: rr_offset_I,  # rr I
                4: rr_offset_Q,  # rr Q
            },
            "digital_outputs": {},
            "analog_inputs": {
                1: {"offset": pars['analog_input_offsets'][0], "gain_db": 3},  # rr I
                2: {"offset": pars['analog_input_offsets'][1], "gain_db": 3}  # rr Q
            },
        },
    },

    "elements": {
        "qubit": {
            "mixInputs": {
                "I": ("con1", 1),
                "Q": ("con1", 2),
                "lo_frequency": qb_LO,
                "mixer": "qubit",
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
                "mixer": "rr",
            },
            "intermediate_frequency": rr_IF,
            "outputs": {
                "out1": ("con1", 1),
                "out2": ("con1", 2),
            },
            "time_of_flight": 264, # should be multiple of 4 (at least 24)
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
            "length": 2*pars['pi_half_len'],
            "waveforms": {
                "I": "pi_wf_i1",
                "Q": "pi_wf_q1",
            },
        },
        "Xpi_pulse": {
            "operation": "control",
            "length": 2*pars['pi_half_len'],
            "waveforms": {
                "I": "pi_wf_i1",
                "Q": "pi_wf_q1",
            },
        },
        "Ypi_pulse": {
            "operation": "control",
            "length": 2*pars['pi_half_len'],
            "waveforms": {
                "I": "pi_wf_q1",
                "Q": "pi_wf_i1",
            },
        },
        "Xpi_half_pulse": {
            "operation": "control",
            "length": pars['pi_half_len'],
            "waveforms": {
                "I": "pi_half_wf_i1",
                "Q": "pi_half_wf_q1",
            },
        },
        "Ypi_half_pulse": {
            "operation": "control",
            "length": pars['pi_half_len'],
            "waveforms": {
                "I": "pi_half_wf_q1",
                "Q": "pi_half_wf_i1",
            },
        },
        "gaussian_pulse": {
            "operation": "control",
            "length": pars['gauss_len'],
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
            "length": pars['pi_half_len'],
            "waveforms": {
                "I": "pi_half_wf_i1",
                "Q": "pi_half_wf_q1",
            },
        },
        "ro_pulse1": {
            "operation": "measurement",
            "length": pars['rr_pulse_len_in_clk'] *4, # in ns (needs to be multiple of 4)
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
        "const_wf": {"type": "constant", "sample": pars['amp_q']},
        "const_wf_rr": {"type": "constant", "sample": pars['amp_r']},
        "gaussian_wf": {"type": "arbitrary", "samples": [float(arg) for arg in pars['gauss_amp'] * gaussian(pars['gauss_len'], pars['gauss_len']/5)]},
        "gaussian_4ns_wf": {"type": "arbitrary", "samples": gauss_wf_4ns},
        "ro_wf1": {"type": "constant", "sample": pars['amp_r']},
        "pi_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_amp'] * gaussian(2*pars['pi_half_len'], 2*pars['pi_half_len']/5)]},
        "pi_wf_q1": {"type": "constant", "sample": 0.0},
        "pi_half_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_amp'] * gaussian(pars['pi_half_len'], pars['pi_half_len']/5)]},
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
        "qubit": [{"intermediate_frequency": qb_IF, "lo_frequency": pars['qb_LO'], "correction": IQ_imbalance(*pars['qubit_mixer_imbalance'])}],
        "rr": [{"intermediate_frequency": rr_IF, "lo_frequency": pars['rr_LO'], "correction": IQ_imbalance(*pars['rr_mixer_imbalance'])}],
    }
}

# qmm = QuantumMachinesManager()
# qm = qmm.open_qm(config)
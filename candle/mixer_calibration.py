# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 11:35:33 2022

Description: Functions used to measure power at a given frequency using spectrum
analyzer and perform mixer calibration. The spectrum analyzer should be initialized
using the <instrument_init.py> script.

@author: Evangelos Vlachos <evlachos@usc.edu>
"""

from qm.QuantumMachinesManager import QuantumMachinesManager
import time
import numpy as np
import scipy.optimize as opti
from config import *
from VISAdrivers.sa_api import *
import matplotlib.pyplot as plt
from qm.qua import *
from tqdm import tqdm
from measurement_functions import play_pulses
from meas_utilities import IQ_imbalance
from plot_functions import plot_mixer_opt
import json
from instrument_init import set_qb_LO,set_rr_LO,set_attenuator,get_attenuation

'''Qubit mixer calibration'''

mixer1 = 'qubit'
# get LO leakage power
qb_lo_leakage = get_power(sa,freq=qb_LO,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get image leakage power
qb_im_leakage = get_power(sa, freq=qb_LO-qb_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get qubit drive power at (almost) top of fridge
qb_on_power = get_power(sa, freq=qb_LO+qb_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
opt_mixer(sa, cal='LO', freq_span = 1e6,  mode='coarse',element=mixer1)
# do a finer sweep in desired range
opt_mixer(sa, cal='LO', freq_span = 1e6, reference = -50, mode='fine',element=mixer1)
# do a coarse sweep to minimize sideband
opt_mixer(sa, cal='SB', freq_span = 1e6,mode='coarse',element=mixer1)
# do a finer sweep in desired range
opt_mixer(sa, cal='SB', freq_span = 1e6, reference = -50, mode='fine',element=mixer1)

'''Do IQ imbalance calibration'''

mixer2 = 'rr'
# get LO leakage power
rr_lo_leakage = get_power(sa, freq=rr_LO,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get image leakage power
rr_im_leakage = get_power(sa, freq=rr_LO-rr_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get qubit drive power at (almost) top of fridge
rr_on_power = get_power(sa, freq=rr_LO+rr_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# do a coarse sweep to minimize LO leakage
opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',element=mixer2)
# do a finer sweep in desired range
opt_mixer( sa, cal='LO',  freq_span = 1e6, reference = -50, mode='fine',element=mixer2)
# do a coarse sweep to minimize sideband
opt_mixer( sa, cal='SB', freq_span = 1e6,mode='coarse',element=mixer2)
# do a finer sweep in desired range
opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = -50, element=mixer2)

sa_close_device(sa)

#%% config_sa
def config_sa(sa,freq,span=5e6,reference=-30):
    """
    Prepares spectrum analyzer for measurement

    Parameters
    ----------
    sa :
        Handle for spectrum analyzer
    freq : float
        Center frequency of span.
    span : float, optional
        DESCRIPTION. The default is 5e6.
    reference : float, optional
        Upper power threshold of SA in dBm. The default is -30.

    Returns
    -------
    None.

    """

    sa_config_level(sa, reference) # sets sensitivity
    sa_config_center_span(sa, freq, span) # sets center frequency
    sa_initiate(sa, SA_SWEEPING, 0)
    query = sa_query_sweep_info(sa)
    sweep_length = query["sweep_length"]
    start_freq = query["start_freq"]
    bin_size = query["bin_size"]
    freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)

#%% get_power
def get_power(sa,freq,reference=-100,span=1e6, config=False, plot=False):
    """
    Configures SA (optional) and measures power at specified frequency

    Parameters
    ----------
    sa : ???
        spectrum analyzer.
    freq : TYPE
        DESCRIPTION.
    reference : TYPE, optional
        DESCRIPTION. The default is -100.
    span : TYPE, optional
        DESCRIPTION. The default is 5e6.
    config : boolean, optional
        whether to reconfigure the SA or not. Set to false when calibrating mixers. The default is False.
    plot : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    freqs : TYPE
        DESCRIPTION.
    power : TYPE
        DESCRIPTION.

    """
    # skips configuring the spectrum analyzer. Used only when optimizing mixer
    if config:
        play_pulses()
        sa_config_level(sa, reference) # sets sensitivity
        sa_config_center_span(sa, freq, span) # sets center frequency
        sa_initiate(sa, SA_SWEEPING, 0)
        query = sa_query_sweep_info(sa)
        sweep_length = query["sweep_length"]
        start_freq = query["start_freq"]
        bin_size = query["bin_size"]
        freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)

    # measure
    signal = sa_get_sweep_64f(sa)['max']
    power = round(np.max(signal),1)

    if plot:
        plt.plot(freqs*1e-6, signal,'-')
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Power [dBm]')
        plt.show()
        print(f'{power} dBm at {freq/1e9} GHz')

    return power

#%% opt_mixer
def opt_mixer(sa,cal,mode,element,freq_span=1e6,reference = -30, plot=True):
    """
    Minimizes leakage at LO ('lo' option) or at image sideband ('sb' option) by sweeping the relevant parameters

    Args:
        sa (): spectrum analyzer handle.
        cal (str): Whether to minimize LO leakage or image sideband. The default is 'lo'.
        mode (str): Coarse of fine stepsize.
        element (str):       Which element to optimize (qubit (qubit) or readout (rr)).
        pars (dict):        dictionary containing experimental values like mixer calibration offsets.
        freq (float):      Frequency at which we want to minimize power.
        reference (float): Threshold of spectrum analyzer.
        plot (TYPE, optional): DESCRIPTION. Defaults to False.

    Returns:
        values (float array): optimal calibration values.
        argmin (TYPE): DESCRIPTION.

    """
    print('Opening parameter JSON file')
    with open('pars.json', 'r') as openfile:
        pars = json.load(openfile)


    qm = play_pulses() # used to generate siebands
    # gets frequency values from config file
    freqLO = config['mixers'][f'{element}'][0]['lo_frequency']
    freqIF = config['mixers'][f'{element}'][0]['intermediate_frequency']

    if cal == 'LO':
        freq = freqLO
        par1 = pars[f'{element}_mixer_offsets'][0]
        par2 = pars[f'{element}_mixer_offsets'][1]
        print(f'LO at {round(freqLO*1e-9,5)} GHz\nCurrent I_offset = {round(par1*1e3,1)} mV, Current Q_offset = {round(par2*1e3,1)} mV')
    elif cal == 'SB':
        freq = freqLO - freqIF
        par1 = pars[f'{element}_mixer_imbalance'][0]
        par2 = pars[f'{element}_mixer_imbalance'][1]
        print(f'Sideband at {round((freqLO-freqIF)*1e-9,5)} GHz\nCurrent gain = {round(par1,4)}, Current phase = {round(par2,4)}')

    if element == 'rr':
        atten = get_attenuation()
        set_attenuator(attenuation=0)

    config_sa(sa,freq,span=freq_span,reference=reference) # setup spectrum analyzer for measurement

    # initialize sweep parameters
    if cal == 'LO':
        if mode == 'coarse':
            span=20e-3
            step=2e-3
        elif mode == 'fine':
            span=2e-3
            step=0.1e-3
    elif cal == 'SB':
        if mode == 'coarse':
            span=150e-3
            step=12.5e-3
        elif mode == 'fine':
            span=75e-3
            step=5e-3

    par1_arr = np.arange(par1-span/2, par1+span/2, step)
    par2_arr = np.arange(par2-span/2, par2+span/2, step)
    L1 = len(par1_arr)
    L2 = len(par2_arr)
    power_data = np.zeros((L1,L2))

    # sweep parameters and get power at every point
    with tqdm(total = L1*L2) as progress_bar:
        for i, par1 in enumerate((par1_arr)):
            for j, par2 in enumerate((par2_arr)):
                if cal == 'LO':
                    qm.set_output_dc_offset_by_element(element, "I", par1)
                    qm.set_output_dc_offset_by_element(element, "Q", par2)
                elif cal == 'SB':
                    qm.set_mixer_correction(element,int(freqIF), int(freqLO), IQ_imbalance(par1, par2))
                power_data[i,j] = get_power(sa, freq,span=freq_span)
                progress_bar.update(1)

    argmin = np.unravel_index(power_data.argmin(), power_data.shape)
    # print(argmin)
    # set the parameters to the optimal values and modify the JSON dictionary
    if cal == 'LO':
        opt_I = par1_arr[argmin[0]]
        opt_Q = par2_arr[argmin[1]]
        qm.set_output_dc_offset_by_element(element, "I", opt_I)
        qm.set_output_dc_offset_by_element(element, "Q", opt_Q)
        pars[f'{element}_mixer_offsets'] = [opt_I,opt_Q]
        with open("pars.json", "w") as outfile:
            json.dump(pars, outfile)
        print(f'optimal I_offset = {round(opt_I*1e3,1)} mV, optimal Q_offset = {round(1e3*opt_Q,1)} mV')
    elif cal == 'SB':
        qm.set_mixer_correction(element,int(freqIF), int(freqLO), IQ_imbalance(par1_arr[argmin[0]],par2_arr[argmin[1]]))
        pars[f'{element}_mixer_imbalance'] = [par1_arr[argmin[0]],par2_arr[argmin[1]]]
        with open("pars.json", "w") as outfile:
            json.dump(pars, outfile)
        print(f"optimal gain = {round(par1_arr[argmin[0]],4)}, optimal phi = {round(par2_arr[argmin[1]],4)}")

    if element == 'rr':
        set_attenuator(atten)

    print(f'Power: {np.amin(power_data)} dBm at {freq/1e9} GHz')
    if plot:
        plot_mixer_opt(par1_arr, par2_arr, power_data,cal=cal,element=element,fc=freq)

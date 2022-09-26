# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 11:35:33 2022

@author: haimeng zhang <haimeng@usc.edu>
"""
from VISAdrivers.sa_api import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from config import *
from qm.qua import *
from tqdm import tqdm

default_sweepBW = 300
default_span = 50e3
peak_dist= default_span/2
low_ref = -130 # baseline of noise level of SA
high_ref = -120
on_power = 0

class mixer_calibration:
    def __init__(self, qm, mixer):
        self.sa = sa_open_device()["handle"]
        self.reference = {'off': high_ref,
                          'on': on_power} # reference level of sa
        if mixer == 'qubit':
            self.element = "qubit"
            self.LO_freq = qb_LO
            self.IF_freq = qb_IF
            pulse = "const"
            self.mixer = "mixer_q1"
        elif mixer == 'readout':  
            self.element = "rr"
            self.LO_freq = rr_LO
            self.IF_freq = rr_IF
            pulse = "readout"
            self.mixer = "mixer_rl1"
            
        with program() as mixer_cal:
            with infinite_loop_():
                play(pulse, self.element, duration=100)
        
        self.qm = qm
        self.job = qm.execute(mixer_cal)
    
    def read_off_leakage_power(self,reference=-100,span=default_span, plot=False):
        marker = self.LO_freq
        sa_config_acquisition(device = self.sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
        sa_config_level(self.sa, reference)
        sa_config_gain_atten(self.sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
        sa_config_sweep_coupling(device = self.sa, rbw = default_sweepBW, vbw = default_sweepBW, reject=0)
        sa_config_center_span(self.sa, marker, span)
        sa_initiate(self.sa, SA_SWEEPING, 0)
        query = sa_query_sweep_info(self.sa)
        sweep_length = query["sweep_length"]
        start_freq = query["start_freq"]
        bin_size = query["bin_size"]
        freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)
        if not freqs[0] <= marker <= freqs[-1]:
            raise ValueError('marker not in frequency range')
        signal = sa_get_sweep_64f(self.sa)
        # offset_marker = abs(freqs-marker)
        try:
            index = self.find_peak(signal['max'], height=0.9*max(signal['max']),distance=span/2)
            power = signal['max'][index[0]]
            bPeak = True
        except IndexError:
            print('peak not found')
            bPeak = False
            pass
        if plot:
            plt.plot(freqs, signal['max'],'-')
            if bPeak:
                plt.plot(freqs[index], signal['max'][index], 'x')
            plt.show()
        return freqs, signal['max']
        
    
    def read_off_image_power(self,reference=-100,span=default_span, plot=False):
        marker = self.LO_freq  - self.IF_freq
        sa_config_acquisition(device = self.sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
        sa_config_level(self.sa, reference)
        sa_config_gain_atten(self.sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
        sa_config_sweep_coupling(device = self.sa, rbw = default_sweepBW, vbw = default_sweepBW, reject=0)
        sa_config_center_span(self.sa, marker, span)
        sa_initiate(self.sa, SA_SWEEPING, 0)
        query = sa_query_sweep_info(self.sa)
        sweep_length = query["sweep_length"]
        start_freq = query["start_freq"]
        bin_size = query["bin_size"]
        freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)
        if not freqs[0] <= marker <= freqs[-1]:
            raise ValueError('marker not in frequency range')
        signal = sa_get_sweep_64f(self.sa)
        # offset_marker = abs(freqs-marker)
        # try:
        #     index = self.find_peak(signal['max'], height=0.95*max(signal['max']),distance=span/2)
        #     power = signal['max'][index[0]]
        #     bPeak = True
        # except IndexError:
        #     print('peak not found')
        #     bPeak = False
        #     pass
        bPeak = True
        index = np.argmax(signal['max'])
        if plot:
            plt.clf()
            plt.plot(freqs, signal['max'],'-')
            if bPeak:
                plt.plot(freqs[index], signal['max'][index], 'x')
            plt.show()
        return freqs, signal['max']
    
    
    def config_sa_sweep(self,center_freq,reference=-100,span=default_span,sweepBW=default_sweepBW):
        sa_initiate(self.sa, SA_SWEEPING, 0)
        sa_config_level(self.sa, reference)
        sa_config_gain_atten(self.sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
        sa_config_sweep_coupling(device = self.sa, rbw = sweepBW, vbw = sweepBW, reject=0)
        sa_config_center_span(self.sa, center_freq, span)
        sa_initiate(self.sa, SA_SWEEPING, 0)
        
    def sweep_and_plot(self,center_freq,reference=-100,span=default_span,sweepBW=default_sweepBW, **kwargs):
        self.config_sa_sweep(center_freq,reference,span,sweepBW)
        query = sa_query_sweep_info(self.sa)
        sweep_length = query["sweep_length"]
        start_freq = query["start_freq"]
        bin_size = query["bin_size"]
        freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)
        if not freqs[0] <= center_freq <= freqs[-1]:
            raise ValueError('marker not in frequency range')
        signal = sa_get_sweep_64f(self.sa)
        if 'ax' in kwargs:
            ax = kwargs.get('ax')
            ax.plot(freqs, signal['max'],'-')
        else:
            plt.plot(freqs, signal['max'],'-')
        try:
            index = self.find_peak(signal['max'], self.reference['off'],distance=peak_dist)
            power = signal['max'][index]
            if 'ax' in kwargs:
                ax = kwargs.get('ax')
                ax.plot(freqs[index], signal['max'][index], 'x')
            else:
                plt.plot(freqs[index], signal['max'][index], 'x')
            return power
        except:
            pass
        plt.show()
        

        
    def find_peak(self, signal, height, distance):
        indices, pdict = find_peaks(signal, height=max(signal),distance=distance)
        if len(indices) > 1:
            indices = sorted(indices, key=lambda x: signal[x],reverse=True)
        return indices[0]
    
    def get_leakage(self, i0, q0):
        self.qm.set_dc_offset_by_qe(self.element, "I", i0)
        self.qm.set_dc_offset_by_qe(self.element, "Q", q0)
        amp_ = self.get_amp()
        return amp_
    
    def brute_force_search_dc_offsets(self, Imin = -0.05, Imax = 0.05, Qmin = -0.05, Qmax = 0.05, num_of_points = 30, plot=False):
     
        dc_offsets_I = np.linspace(Imin, Imax, num_of_points)
        dc_offsets_Q = np.linspace(Qmin, Qmax, num_of_points)
        values = np.zeros((num_of_points, num_of_points))
        
        with tqdm(total = num_of_points**2) as progress_bar:
            
            for i, i_value in enumerate((dc_offsets_I)):
                
                for j, q_value in enumerate((dc_offsets_Q)):
            
                    values[i,j] = self.get_leakage(i_value, q_value)
                    progress_bar.update(1)
                    
        argmin = np.unravel_index(values.argmin(), values.shape)
        print(f'optimal i_offset = {dc_offsets_I[argmin[0]]}, optimal q_offset = {dc_offsets_Q[argmin[1]]}')
        
        if plot:
            plt.clf()
            im = plt.imshow(values, extent=(Qmin, Qmax, Imin, Imax), origin='lower')
            plt.colorbar(im)
            plt.xlabel("Q")
            plt.ylabel("I")
            plt.show()
        
        return values, argmin
    
    def IQ_imbalance_correction(self, g, phi):
        c = np.cos(phi)
        s = np.sin(phi)
        N = 1 / ((1 - g ** 2) * (2 * c ** 2 - 1))
        return [
            float(N * x) for x in [(1 - g) * c, (1 + g) * s, (1 - g) * s, (1 + g) * c]
        ]
    
    def get_image(self, g, p):
        self.qm.set_mixer_correction(self.mixer,int(self.IF_freq), int(self.LO_freq), self.IQ_imbalance_correction(g, p))
        amp_ = self.get_amp()
        return amp_
    
    
    def brute_force_search_imbalance(self, gmin=-0.15, gmax=0.15, pmin=-0.15, pmax=0.15, num_of_points=30, plot=False):
        
        imbalances_g = np.linspace(gmin, gmax, num_of_points)
        imbalances_p = np.linspace(pmin, pmax, num_of_points)
        values = np.zeros((num_of_points, num_of_points))
        
        with tqdm(total = num_of_points**2) as progress_bar:
        
            for i, g_value in enumerate(imbalances_g):
                
                for j, phi_value in enumerate(imbalances_p):
                    values[i,j] = self.get_image(g_value, phi_value)
                    progress_bar.update(1)
                
        argmin = np.unravel_index(values.argmin(), values.shape)
        print(f"optimal gain = {imbalances_g[argmin[0]]}, optimal phi = {imbalances_p[argmin[1]]}")
        if plot:
            plt.clf()
            im = plt.imshow(values, extent=(pmin, pmax, gmin, gmax), origin='lower')
            plt.colorbar(im)
            plt.xlabel("phase")
            plt.ylabel("gain")
            plt.show()
        return values, argmin
    
    def get_amp(self):
        signal = sa_get_sweep_64f(self.sa)
        index = self.find_peak(signal['max'], height=0.9*max(signal['max']),distance=default_span/2)
        return signal['max'][index]
        
        
    
    def close(self):
        self.job.halt()
        sa_close_device(self.sa)
        

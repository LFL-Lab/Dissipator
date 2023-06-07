# -*- coding: utf-8 -*-
"""
Created on Thu May 11 15:07:34 2023

@author: lfl

data analysis code for ring down 
"""
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit

def decay(x,amp,tau,offset):
    return amp*np.exp(-x/tau)+offset

class Experiment:
    
    def __init__(self, dataDir, filename):
        self.dataDir = dataDir
        self.filename = filename
        f = h5py.File(self.dataDir + self.filename, 'r')
        self.data_key_list = list(f.keys())
        self.sweep_list = {}
        self.data = {}
        f.close()
        
        
    def read_sweep_parameters(self):
        f = h5py.File(self.dataDir + self.filename, 'r')
        data_key = list(f.keys())[1]
        freq_list = list(set([k.split('ffl_freq=')[1].split('GHz')[0] for k in f[data_key].keys()]))
        flux_list = list(set([k.split('flux=')[1].split('uA')[0] for k in f[data_key].keys()]))
        freq_list.sort()
        flux_list.sort()
        f.close()
        self.sweep_list['freq_list'] = freq_list
        self.sweep_list['flux_list'] = flux_list
        
    def plot_sweep_data(self, f, key, plot=True):
        flux = key.split('flux=')[1].split('uA')[0]
        fflFreq = key.split('ffl_freq=')[1].split('GHz')[0]
        amp_keys = list(f.keys())
       
        num_of_entries = len(amp_keys)
        if plot:
            fig, axes = plt.subplots(1, 2, figsize=(8,5))
    
        t1 = np.zeros((num_of_entries,))
        t1_err = np.zeros((num_of_entries,))
    
        cmap =['#000000',
                '#e69f00',
                '#56b4e9',
                '#009e73',
                '#cc79a7',
                '#f0e442',
                '#0072b2',
                '#d55e00',
                ]
        colors = sns.color_palette("Set3", num_of_entries)
    
        for i in range(num_of_entries):
            entry0 = amp_keys[i]
            time = np.array(f[entry0]['time'])
            data = np.array(f[entry0]['I']) + 1j * np.array(f[entry0]['Q'])
            mag = abs(data) * 1e3
            tau = 0.2
            amp = abs(mag[0] - mag[-1])
            offset = 0
    
            p0 = [amp,tau,offset]
            fitFunction = decay
            fitted_pars, covar = curve_fit(fitFunction, time, mag,p0=p0,method='trf',xtol=1e-12,maxfev=40e3)
            error = np.sqrt(abs(np.diag(covar)))
            t1[i] = fitted_pars[1]
            t1_err[i] = error[1]
    
            if plot:
                axes[0].plot(time, mag ,'x',color = colors[i])
                if i != 0:
                    label='drive amp='+entry0.split(' = ')[1] + ', '
                else:
                    label = 'drive off, '
    
                axes[0].plot(time, decay(time,fitted_pars[0], fitted_pars[1], fitted_pars[2]),'-',color =colors[i],label=label+r'$\tau$'+f'={fitted_pars[1]*1e3:.0f}'+r'$\pm$'+f'{t1_err[i]*1e3:.0f}'+r'$\mu$s')
                axes[0].set_xlabel('delay time ('+r'$\mu$'+'s)')
                axes[0].set_ylabel('voltage (mV)')
                axes[0].legend(fontsize=7,bbox_to_anchor=(1.3, -0.22),ncol=2)
            axes[0].set_xlim(0, 1)
                
        if plot:
            amp_list = [round(0.0 + 0.1* n,1) for n in range(len(amp_keys)-1)]
            self.sweep_list['amp_list'] = amp_list
            axes[-1].set_xlabel('drive amplitude scaling')
            axes[-1].errorbar(amp_list, t1[1:], yerr=t1_err[1:],color=colors[0], label=f'{flux}uA, {fflFreq}GHz')
            axes[-1].legend()
            fig.tight_layout()
            fig.show()
        
        return t1, t1_err
    
    def analyze_data(self, data_key, plot=True):
        f = h5py.File(self.dataDir + self.filename, 'r')
        num_of_entries = 12
        num_of_fluxes = len(self.sweep_list['flux_list'])
        num_of_freqs = len(self.sweep_list['freq_list'])
        data = np.zeros((num_of_freqs, num_of_entries, num_of_fluxes))
        data_err = np.zeros((num_of_freqs, num_of_entries, num_of_fluxes))
        
        
        for i,flux in enumerate(self.sweep_list['flux_list']):
            for j,freq in enumerate(self.sweep_list['freq_list']):
                k = [s for s in f[data_key].keys() if f'flux={flux}uA_ffl_freq={freq}GHz_' in s][0]
                group = f[data_key][k]
                t1, t1_err = self.plot_sweep_data(group, k, plot=plot)
                data[j,:,i] = t1
                data_err[j,:,i] = t1_err
        self.data['t1'] = np.stack((data, data_err), axis=-1)
        f.close()
        
    def plot_data(self):
        self.sweep_list['freq_list'] = np.array([float(f) for f in self.sweep_list['freq_list']])
        fig,ax = plt.subplots(1,1, figsize=(6,8))
        stop_freq = 5
        stop_index = np.argmin(abs(np.array(self.sweep_list['freq_list'])-stop_freq))
        for i,flux in enumerate(self.sweep_list['flux_list']):
            flux_data = self.data['t1'][:,:,i,0]
            plt.rcParams.update({'font.size': 12})
            plt.imshow(flux_data[:stop_index,1:], cmap='RdBu', interpolation='nearest',origin='lower', 
                       extent=[self.sweep_list['amp_list'][0],self.sweep_list['amp_list'][-1], 
                               self.sweep_list['freq_list'][0], self.sweep_list['freq_list'][stop_index]],
                       vmin=0, vmax=0.5)
            plt.title(f'flux={flux}uA')
            plt.colorbar()
            plt.xlabel('drive amplitude scaling')
            plt.ylabel('drive frequency (GHz)')
            plt.show()
        return fig
#%% plotting data
res_flux = Experiment(dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss09_6024\\20230518\\fluxSweep\\',
                filename = 'fluxSweep_diss09_6024_start=-200uA_stop=200uA_5.h5')
res_flux.plot_data()

#%% example
t1 = Experiment(dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss09_6024\\20230515\\ringdown\\',
                filename = 'ringdown_sweepFLux_sweepFFLfreq_ffl_IF=0.0MHz_amp_ffl=0.3_DA=24dB_fDA=20dB_navg=2000_2.h5')
t1.read_sweep_parameters()
t1.analyze_data(data_key=t1.data_key_list[1])
fig = t1.plot_data()

# -*- coding: utf-8 -*-
"""
Created on Sat May  6 00:14:24 2023

@author: lfl

resonator spectroscopy as a function of flux and calibrate FFL flux modulation amplitude
"""
import timeit
from dissipator import *
from sequence import *
import instrument_init as inst
import plot_functions as pf
import h5py
from datetime import datetime
import os
from instrument_init import init_sa, init_sa_by_serial_number


def acquire_rr_spec_background(qb, base_flux = 56e-6, 
                               IF_min=40e5, 
                               IF_max=54e6,
                               n_avg=1000,):
    inst.set_flux_bias(base_flux)
    I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=IF_min,IF_max=IF_max,df=0.1e6,n_avg=n_avg,savedata=True,fit=False, flux=base_flux)
    dataDict = {'metadata': {'flux': base_flux,
                             'report': str(job.execution_report())},
                'freqs': freqs,
                'I': I,
                'Q': Q,}
    return dataDict

def main():
    device = 'diss09_6024'
    qb = dissipator(device, device_name=device)
    '''Update important parameters'''
    qb.update_value('rr_freq', 6.025539e9)
    qb.update_value('rr_LO', 5.975e9)
    qb.update_value('rr_IF', qb.pars['rr_freq'] -qb.pars['rr_LO'] )
    qb.update_value('rr_atten', 24)
    n_avg = 4000
    # save data
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\calibration'
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    filename = f'fflCalibration_DA={qb.pars["rr_atten"]}dB_navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)
    # move flux, compute how much resonance moved
    base_flux = 51e-6
    start_flux = -80e-6
    stop_flux = 5e-6
    step_size = 5e-6
    flux_list = np.arange(start_flux, stop_flux+step_size/2, step_size)
    coil_fc_dict = {}
    ffl_fc_dict = {}
    
    inst.turn_on_ffl_source_meter()
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_coil = hf.create_group(f'rrSpecCoil_{timestamp}')
        g_ffl = hf.create_group(f'rrSpecFFL_{timestamp}')
        base_flux = 51e-6
        dataDict = acquire_rr_spec_background(qb, base_flux=base_flux)
        save_datadict_to_fgroup(g_coil, f'base flux = {base_flux*1e6:.0f} uA', dataDict)
        I0 = dataDict['I']
        Q0 = dataDict['Q']
        freqs = dataDict['freqs']
        z = np.polyfit(freqs, np.abs(I0 + 1j*Q0), 3)
        p = np.poly1d(z)
        for flux in tqdm(flux_list):
            inst.set_flux_bias(flux)
            I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=40e6,IF_max=54e6,df=0.1e6,n_avg=n_avg,savedata=True)
            fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q) - p(freqs))
            # plot 
            pf.spec_plot(freqs,I,Q,attenuation=qb.pars['rr_atten'],df=0.1e6,element='resonator',fwhm=fwhm,fc=fc, flux=flux)
            coil_fc_dict[f'coil flux = {flux*1e6:.0f} uA']= fc
            dataDict = {
                'I': I,
                'Q': Q,
                'freqs': freqs,
                'metadata': {'flux': flux,
                             'rr_freq': fc},
                }
            save_datadict_to_fgroup(g_coil, f'coil flux = {flux*1e6:.0f} uA', dataDict)
        print(f'coil flux from {start_flux*1e6:.0f} to {stop_flux*1e6:.0f} uA, resonance has moved from {coil_fc_dict[f"coil flux = {start_flux*1e6:.0f} uA"]/1e9:.6f} to {coil_fc_dict[f"coil flux = {stop_flux*1e6:.0f} uA"]/1e9:.6f} GHz')


    
    # flux_list_reverse = np.flip(flux_list)
    start_flux = 0e-6
    stop_flux = 1500e-6
    step_size = 30e-6
    flux_list_reverse = np.arange(start_flux, stop_flux+step_size/2, step_size)
    
    # use FFL to move resonance back
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','a') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_ffl = hf.create_group(f'rrSpecFFL_{timestamp}')
        
        for flux in tqdm(flux_list_reverse):
            inst.set_ffl_bias(flux, upper_bound=flux + 5e-6)
            I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=40e6,IF_max=54e6,df=0.1e6,n_avg=n_avg,savedata=True,fit=False)
            fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q) - p(freqs))
            # plot 
            pf.spec_plot(freqs,I,Q,attenuation=qb.pars['rr_atten'],df=0.1e6,element='resonator',fwhm=fwhm,fc=fc, flux=flux)
            ffl_fc_dict[f'ffl flux = {flux*1e6:.0f} uA'] = fc
            dataDict = {
                'I': I,
                'Q': Q,
                'freqs': freqs,
                'metadata': {'flux': flux,
                             'rr_freq': fc},
                }
            save_datadict_to_fgroup(g_ffl, f'ffl flux = {flux*1e6:.0f} uA', dataDict)
        print(f'ffl flux from {flux_list_reverse[0]*1e6:.0f} to {flux_list_reverse[-1]*1e6:.0f} uA, resonance has moved from {ffl_fc_dict[f"ffl flux = {start_flux*1e6:.0f} uA"]/1e9:.6f} to {ffl_fc_dict[f"ffl flux = {stop_flux*1e6:.0f} uA"]/1e9:.6f} GHz')
    
    for k,v in coil_fc_dict.items():
        print(f'{k}:{v}')
    for k,v in ffl_fc_dict.items():
        print(f'{k}:{v}')
if __name__ == "__main__":
    main()
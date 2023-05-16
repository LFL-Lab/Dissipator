# -*- coding: utf-8 -*-
"""
Created on Sun May 14 15:43:09 2023

@author: lfl

collect flux sweep data
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

qb = dissipator('diss09_6024', device_name='diss09_6024')
qb.update_value('diss_freq', 10.1e9)
qb.update_value('rr_LO', 5.975e9)
qb.update_value('rr_atten', 24)
inst.turn_off_ffl_drive()
inst.set_rr_LO(qb.pars['rr_LO']) # turn on
n_avg = 4000
bOptimizeFFLMixer = True
bOptimizeRRMixer = False
bCalibrateRo = True

start_flux = -150e-6
stop_flux = 150e-6
step_size = 20e-6
flux_list = np.arange(start_flux, stop_flux+step_size/2, step_size)
# resonator spectroscopy span is 25 MHz
start_freq = 6.0125e9
stop_freq = 6.0375e9
IF_min = start_freq - qb.pars['rr_LO']
IF_max = stop_freq - qb.pars['rr_LO']
today = datetime.today()
sDate =  today.strftime("%Y%m%d")
saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\spectroscopy\\rr_spec'
if not os.path.exists(saveDir):
    Path(saveDir).mkdir(parents=True, exist_ok=True)
filename = f'rrSpecFluxSweep_DA={qb.pars["rr_atten"]}dB_navg={n_avg}'
index = get_index_for_filename(saveDir, filename)

with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
    now = datetime.now()
    timestamp = now.strftime("%H:%M:%S")
    g_rr = hf.create_group(f'rrSpecCoil_{timestamp}')
    for flux in tqdm(flux_list):
        inst.set_flux_bias(flux,lower_bound=-155e-6,upper_bound=155e-6)
        I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=36e6,IF_max=65e6,df=0.1e6,n_avg=n_avg,savedata=True)
        dataDict = {
            'I': I,
            'Q': Q,
            'freqs': freqs,
            'metadata': {'flux': flux,
                        },
            }
        save_datadict_to_fgroup(g_rr, f'flux = {flux*1e6:.0f} uA', dataDict)
 
inst.set_flux_bias(0.0)


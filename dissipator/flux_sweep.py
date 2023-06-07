

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
qb = dissipator('diss09b_611',device_name='diss09b_611')
#qb.update_value('diss_freq', 10.1e9)
#qb.update_value('rr_LO', 5.975e9)
qb.update_value('rr_atten', 30)
inst.turn_off_ffl_drive()
inst.set_rr_LO(qb.pars['rr_LO']) # turn on
n_avg = 5000
bOptimizeFFLMixer = False
bOptimizeRRMixer = False
bCalibrateRo = False
stop_flux = 150e-6
start_flux = -150e-6
step_size = 10e-6
flux_list = np.arange(start_flux, stop_flux+step_size/2, step_size)
# resonator spectroscopy span is 25 MHz
start_freq = 6.109e9
stop_freq = 6.113e9
IF_min = start_freq - qb.pars['rr_LO']
IF_max = stop_freq - qb.pars['rr_LO']
today = datetime.today()
sDate =  today.strftime("%Y%m%d")
saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\spectroscopy\\rr_spec'
if not os.path.exists(saveDir):
    Path(saveDir).mkdir(parents=True, exist_ok=True)
filename = f'rrSpecFluxSweep_DA={qb.pars["rr_atten"]}dB_navg={n_avg}'
index = get_index_for_filename(saveDir, filename)
df=0.025e6
span=IF_max-IF_min
I = np.zeros((len(flux_list),int(span/df)+1))
Q = np.zeros((len(flux_list),int(span/df)+1))
mag=np.zeros((len(flux_list),int(span/df)+1))
j=0
for flux in tqdm(flux_list):
    inst.set_flux_bias(flux, lower_bound=-300e-6, upper_bound=1e-3)
    dataI, dataQ, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=40,IF_min=IF_min,IF_max=IF_max,df=df,n_avg=10000,savedata=True)
    I[j,:] = dataI
    Q[j,:] = dataQ
    mag[j,:] = np.abs(dataI+1j*dataQ)
    dataDict = {
        'I': I,
        'Q': Q,
        'freqs': freqs,
        'metadata': {'flux': flux,
                        },
        }
    #save_datadict_to_fgroup(g_rr, f'flux = {flux*1e6:.0f} uA', dataDict)
    j += 1
pf.heatplot(xdata=np.around(freqs*1e-9,4),ydata=flux_list*1e6,data=pf.Volt2dBm(mag),xlabel='Frequency (GHz)',ylabel='FLux',cbar_label='Magnitude (dBm)')
inst.set_flux_bias(41e-6,lower_bound=-300e-6,upper_bound=1e-3)

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 17:40:15 2023

@author: lfl
"""

"""
resonator spectroscopy as a function of flux bias
"""

from dissipator import *
import instrument_init as inst
device = 'diss09_5578'
qb = dissipator('diss09', device_name=device)
qb.update_value('rr_LO', 5.5e9)
inst.set_rr_LO(qb.pars['rr_LO'])
start_flux = -100e-6
stop_flux = 800e-6
step = 5e-6
flux_list = np.round(np.arange(start_flux, stop_flux + step/2, step),7)
center = 5.58e9
span = 25e6
num_of_points=1601
start = int(center-span/2)
stop = int(center+span/2)
df = int(span/num_of_points-1)
IF_min = start - qb.pars['rr_LO']
IF_max = stop - qb.pars['rr_LO']
if IF_min <0 or IF_min > 400e6:
    raise ValueError('IF frequency out of range')
freqs = np.arange(start, stop+df/2, df)
Idata = np.zeros((len(flux_list), len(freqs)))
Qdata = np.zeros((len(flux_list), len(freqs)))


today = datetime.today()
sDate =  today.strftime("%Y%m%d")
saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\fluxSweep'
if not os.path.exists(saveDir):
    Path(saveDir).mkdir(parents=True, exist_ok=True)
now = datetime.now()
timestamp = now.strftime("%H:%M:%S")
filename = f'fluxSweep_{device}_start={round(start_flux*1e6)}uA_stop={round(stop_flux*1e6)}uA'
index = get_index_for_filename(saveDir, filename)
with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
    for i,flux in enumerate(tqdm(flux_list)):
        inst.set_flux_bias(flux)
        I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],
                                             atten=qb.pars['rr_atten'],
                                             IF_min=IF_min,IF_max=IF_max,df=df,n_avg=1000,savedata=True)
        dataDict = {'freqs': freqs,
                    'I': I,
                    'Q': Q}
        save_datadict_to_fgroup(hf, f'flux = {flux*1e6} uA', dataDict)
        Idata[i] = I
        Qdata[i] = Q

inst.set_flux_bias(0.0)
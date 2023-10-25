# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 22:55:14 2023

@author: lfl

sweep rr_atten, amp_q_scaling, amp_r_scaling to resolve photon number splitting
"""

from qubit import *
import instrument_init as inst
import h5py
ref_H = 20
ref_L = -45

device = 'diss08_07A'
today = datetime.today()
sDate =  today.strftime("%Y%m%d")
saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\photonResSpec'


qb = qubit('logical')
optimize_mixer = False
qubit_LO = 4.5e9
qb.update_value('qubit_LO', value = qubit_LO)
inst.set_qb_LO(qubit_LO)
rr_atten = 35 # 

if optimize_mixer:
    '''Qubit mixer calibration
    Get leakage power '''
    qb.play_pulses()
    qb_lo_leakage = qb.get_power(sa, freq=qb.pars['qubit_LO'],reference=ref_H,config=True,plot=True)
    qb_im_leakage = qb.get_power(sa, freq=qb.pars['qubit_LO']-qb.pars['qubit_IF'],reference=ref_H,config=True,plot=True)
    qb_on_power = qb.get_power(sa, freq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
    
    '''Optimize mixer'''
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'qubit')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'qubit')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine', element = 'qubit')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'qubit')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'qubit')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element = 'qubit')
    
    '''Readout mixer calibration'''
    # set DA to 0 dB attenuation
    set_attenuator(0)
    get_attenuation()
    rr_lo_leakage = qb.get_power(sa, freq=qb.pars['rr_LO'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
    rr_im_leakage = qb.get_power(sa, freq=qb.pars['rr_LO']-qb.pars['rr_IF'],span = 1e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
    rr_on_power = qb.get_power(sa, freq=qb.pars['rr_LO']+qb.pars['rr_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
    
    # do a coarse sweep to minimize LO leakage
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element='rr')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
    qb.opt_mixer( sa, cal='LO',  freq_span = 1e6, reference = ref_L, mode='fine',element='rr')
    qb.opt_mixer( sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element='rr')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element='rr')

inst.set_attenuator(qb.pars['rr_atten'])

if not os.path.exists(saveDir):
    Path(saveDir).mkdir(parents=True, exist_ok=True)
rr_atten = 33
qb.update_value('rr_atten', rr_atten)
n_avg = 3000
filename = f'photonNumResSpec_DA={rr_atten}dB_navg={n_avg}'
iteration = get_index_for_filename(saveDir, filename)
with h5py.File(f'{saveDir}\\{filename}_{iteration}.h5','w') as hf:
    now = datetime.now()
    timestamp = now.strftime("%H:%M:%S")
    g_on = hf.create_group(f'rr_atten = {rr_atten}')
    amps = [0.01,0.1,0.35,0.65,0.95]
    for amp in amps:
        qb.qubit_spec(f_LO=4.1e9,
                      amp_q_scaling=amp,
                      IF_min=0,
                      IF_max=200e6,
                      df=200e3,
                      n_avg=3000,
                      on_off=True,
                      showprogress=True,
                      savedata=True,
                      check_mixers=False,
                      )
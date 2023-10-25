# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:42:53 2022

@author: lfl
"""
sa = sa
# from mixer_calibration import get_power,opt_mixer
from dissipator import dissipator
from qubit import qubit
import instrument_init as inst
import plot_functions as pf
import numpy as np
from ringdown_ffl import measure_ringdown_drive_off, measure_ringdown_drive_on
from fitTools.Resonator import Resonator


from t1_ffl import measure_t1_ffl_off, measure_t1_w_ffl, resonator_spec_wffl
from g_e_f_T1 import measure_leakage_w_ffl
ref_H = -10
ref_L = -50
'''Initialize qubit class'''
qb = dissipator('diss08_11a',device_name='diss08_11a')
sa =  inst.init_sa()
sa_close_device(sa)

'''Update important parameters'''
qb.update_value('rr_freq', qb.pars['rr_IF'] +qb.pars['rr_LO'])
qb.update_value('rr_LO', 5.39e9)
qb.update_value('rr_IF', qb.pars['rr_freq'] -qb.pars['rr_LO'] )
qb.update_value('rr_atten', 28)
qb.update_value('diss_freq', 9.5e9)
qb.update_value('qubit_freq', 3.36756e9+0.03e6+0.16e6)
qb.update_value('qubit_LO', 3.1e9)
qb.update_value('qubit_IF',qb.pars['qubit_freq'] - qb.pars['qubit_LO'])
qb.update_value('fflqc_LO',4.79e9)
qb.update_value('fflqc_IF',0)
qb.update_value('fflqc_freq',qb.pars['fflqc_freq'] - qb.pars['fflqc_IF'])

#qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['qubit_freq'])
#qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['rr_freq'])
qb.update_value('ffl_freq',2.955e9)
qb.update_value('ffl_LO',2.955e9-300e6)
qb.update_value('ffl_IF',qb.pars['ffl_freq'] - qb.pars['ffl_LO'] )
#qb.update_value('ffl_IF',0. )
qb.update_value('ffl_freq',6.135e9)
qb.update_value('ffl_LO',5.9e9)
qb.update_value('ffl_IF',qb.pars['ffl_freq'] - qb.pars['ffl_LO'] )





inst.set_rr_LO(qb.pars['rr_LO'])
'''Qubit mixer calibration 
Get leakage power '''
inst.set_qb_LO(qb.pars['qubit_LO'])
qm = qb.play_pulses(element='qubit',amp_scale=1)


qb_lo_leakage = qb.get_power(sa, freq=qb.pars['qubit_LO'],reference=-20,config=True,plot=True)
qb_im_leakage = qb.get_power(sa, freq=qb.pars['qubit_LO']-qb.pars['qubit_IF'],reference=ref_H,config=True,plot=True)
qb_on_power = qb.get_power(sa, freq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power

'''Optimize mixer''' 
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, amp_q=1,mode='coarse', element = 'qubit')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, amp_q=1,mode='intermediate', element = 'qubit')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, amp_q=1,mode='fine', element = 'qubit')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, amp_q=1,mode='coarse', element = 'qubit')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, amp_q=1,mode='intermediate', element = 'qubit')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, amp_q=1,mode='fine',element = 'qubit')

qb.update_value('qubit12_mixer_offsets',qb.pars['qubit_mixer_offsets'])



'''Readout mixer calibration'''
# set DA to 0 dB attenuation
inst.set_rr_LO(qb.pars['rr_LO'])
qm = qb.play_pulses(element='rr')
inst.set_attenuator(0)
inst.get_attenuation()
rr_lo_leakage = qb.get_power(sa, freq=qb.pars['rr_LO'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
rr_im_leakage = qb.get_power(sa, freq=qb.pars['rr_LO']-qb.pars['rr_IF'],span = 1e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
rr_on_power = qb.get_power(sa, freq=qb.pars['rr_LO']+qb.pars['rr_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',amp_q=1,reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='intermediate',amp_q=1,reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L,amp_q=1, mode='fine',element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6,  mode='coarse',amp_q=1, reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='intermediate',amp_q=1,reference = ref_H, element='rr')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', amp_q=1,reference = ref_L, element='rr')


'''FFL mixer calibration'''
ref_H = -10
#ref_L = -30
ref_L = -50
inst.set_ffl_LO(qb.pars['ffl_LO'])
qm = qb.play_pulses(element='ffl', amp_scale=1, switch='on')
 # turn on

qb_lo_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO'],reference=ref_H,config=True,plot=True)
qb_im_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO']-qb.pars['ffl_IF'],reference=ref_L,config=True,plot=True)
qb_on_power = qb.get_power(sa, freq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power

'''Optimize FFL mixer'''
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H,amp_q=1, mode='coarse', element = 'ffl', switch='on')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H,amp_q=1, mode='intermediate', element = 'ffl', switch='on')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L,amp_q=1, mode='fine', element = 'ffl', switch='on')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H,amp_q=1, mode='coarse', element = 'ffl')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H,amp_q=1, mode='intermediate', element = 'ffl')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L,amp_q=1, mode='fine',element = 'ffl')


'''FFLqc mixer calibration'''
ref_H = -10
#ref_L = -30
ref_L = -50
inst.set_ffl_LO(qb.pars['fflqc_LO'])
qm = qb.play_pulses(element='fflqc', switch='off')
 # turn on

qb_lo_leakage = qb.get_power(sa, freq=qb.pars['fflqc_LO'],reference=ref_H,config=True,plot=True)
qb_im_leakage = qb.get_power(sa, freq=qb.pars['fflqc_LO']-qb.pars['fflqc_IF'],reference=ref_H,config=True,plot=True)
qb_on_power = qb.get_power(sa, freq=qb.pars['fflqc_LO']+qb.pars['fflqc_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
'''Optimize fflqc mixer'''
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H,amp_q=1, mode='coarse', element = 'fflqc', switch='off')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H,amp_q=1, mode='intermediate', element = 'fflqc',switch='off')
qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L,amp_q=1, mode='fine', element = 'fflqc',switch='off')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H,amp_q=1, mode='coarse', element = 'fflqc')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H,amp_q=1, mode='intermediate', element = 'fflqc')
qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L,amp_q=1, mode='fine',element = 'fflqc')


#%% time of flight
inst.set_attenuator(0)
qb.tof_cal()

qb.write_pars()
#%% resonator spectroscopy

# I,Q,freqs,job = qb.run_scan(df = 25e3, n_avg = 250, element='resonator', chunksize = 90e6, attenuation=35, lo_min = 5.636e9, lo_max = 5.636e9,
#            showprogress=True, res_ringdown_time = int(4e3))
inst.turn_off_ffl_drive()
inst.turn_off_fflqc_drive()
inst.turn_qb_LO_off()
qb.update_value('qubit_LO', 2.05e9)
inst.set_rr_LO(qb.pars['rr_LO'])
inst.set_qb_LO(qb.pars['qubit_LO'])
inst.set_attenuator(qb.pars['rr_atten'])
#qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=34,IF_min=200e6,IF_max=208e6,df=50e3,n_avg=5000,savedata=True)
inst.set_ffl_bias(1.15e-3, step = 10e-6, lower_bound=-5e-3, upper_bound=6e-3)
I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=31,IF_min=198e6,IF_max=203e6,df=10e3,n_avg=15000,savedata=True)
fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
qb.update_value('rr_freq', fc)
qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'] )
res=Resonator(port_type='notch',f_data=freqs, z_data=(I+1j*Q))
res.autofit()
res.show()
I, Q, freqs, job = qb.resonator_spec_wffl(f_LO=qb.pars['rr_LO'],atten=31,IF_min=200e6,IF_max=206e6,df=50e3,n_avg=5000,savedata=True,)

#%% punchout
I, Q, freq_arr, job = qb.punchout(IF_min=200e6,IF_max=205e6,df=50e3, atten_range=[10,32], atten_step=1, n_avg=1000, f_LO=[qb.pars['rr_LO']])
#%% qubit spectroscopy
qb.update_value('rr_atten', 31)
qb.update_value('ffl_atten', 9)
inst.set_ffl_attenuator(qb.pars['ffl_atten'])
inst.set_ffl_bias(0.80e-3, step = 10e-6, lower_bound=-5e-3, upper_bound=6e-3)
I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=31,IF_min=198e6,IF_max=204e6,df=50e3,n_avg=5000,savedata=True)
fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
qb.update_value('rr_freq', fc)
qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'] )
I,Q,freqs,job = qb.run_scan(df = 0.5e6, n_avg =4000, element='fflqb', chunksize = 350e6,  lo_min = 4.7e9, lo_max = 5.2e9, amp_q_scaling = 1, saturation_dur = 25e3, showprogress=True, res_ringdown_time = int(40e3), flux=inst.get_ffl_bias()*1e3, check_mixers=True)

#qb.update_value('qubit_LO', 2.25e9)
inst.set_qb_LO(qb.pars['qubit_LO'])
I,Q,freqs,job = qb.qubit_spec(sa=sa,
                              f_LO=qb.pars['qubit_LO'],
                              amp_q_scaling=0.0013,amp_r_scaling=1,
                              IF_min=2765e5,IF_max=2774e5,df=10e3,
                              n_avg=10000,
                              rr_freq=qb.pars['rr_freq'],atten=qb.pars['rr_atten'],
                              on_off=True,showprogress=True,savedata=True,check_mixers=False, saturation_dur = int(50e3))






#%% Rabi
inst.turn_on_ffl_drive()
qb.update_value('rr_atten', 30)
inst.set_attenuator(qb.pars['rr_atten'])
qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=35e6,IF_max=65e6,df=0.1e6,n_avg=2000,savedata=True)


qb.update_value('pi_amp', 0.1004)
qb.update_value('pi_half_amp',qb.pars['pi_amp']/2)
amps, I, Q, job, period, fits = qb.power_rabi(a_min = 0.01, a_max = 1., da = 0.01,  n_avg = 2000,fit = True, plot = True, detuning = 0e6, check_mixers = False)

#qb.update_value('pi_amp', 0.098)
#qb.update_value('pi_half_len', 32)

# qb.update_value('qubit_freq', 2.2681e9)
# qb.update_value('qubit_LO', 2.2e9)
# inst.set_qb_LO(qb.pars['qubit_LO'])
# qb.update_value('qubit_IF',qb.pars['qubit_freq'] - qb.pars['qubit_LO'] )

#qb.update_value('pi_amp', 0.4)

t_arr, I, Q, job, pi2width = qb.pulse_exp(exp='rabi',n_avg = 2000, tmin = 40, tmax = 240, dt = 4, amp_q_scaling = 1, fit = True, plot = True, detuning = 0e6, check_mixers=False)
qb.update_value('pi_half_len', 40)
qb.update_value('pi_len', 40)
# plot rabi
x_vector = t_arr
y_vector = np.stack((I,Q),axis=0)
fig, axes = plt.subplots(2,2)
axes = axes.flatten()
axes[0].plot(x_vector*1e3, y_vector[0]*1e3, '-o', markersize = 3, c='C0')
axes[1].plot(x_vector*1e3, y_vector[1]*1e3, '-o', markersize = 3, c='C0')
axes[2].plot(x_vector*1e3, np.angle(y_vector[0]+1j*y_vector[1]), '-o', markersize = 3, c='C0')
axes[3].plot(x_vector*1e3, np.abs(y_vector[0]+1j*y_vector[1])*1e3, '-o', markersize = 3, c='C0')
for ax in axes:
    ax.set_ylabel('Voltage (mV)')
    ax.set_xlabel('Pulse Duration (ns)')
qubitDriveFreq = qb.pars['qubit_freq']
textstr = '$\omega_d$ = %.4f GHz\n$T_{\pi/2}$ = %.1f ns\n$\hat{n}$ = %d'%(qubitDriveFreq*1e-9,0,5000)
fig.tight_layout()
plt.gcf().text(1.0, 0.15, textstr, fontsize=14)
plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)
plt.show()
#%% optimize readout


#%% Ramsey & T1
inst.turn_off_ffl_drive()
#inst.turn_on_ffl_drive()
t_arr, I, Q, job, fitted_pars = qb.pulse_exp(exp='T1', check_mixers = False, n_avg=2000,dt=1024,tmin=0,tmax=100e3,detuning=0e3, amp_ffl_scale=0, flux=inst.get_ffl_bias()*1e3)
t_arr, I, Q, job, fitted_pars = qb.pulse_exp(exp='ramsey', check_mixers = False, n_avg=5000,dt=256,tmin=0,tmax=12e3,detuning=0.16e6)


#%% Echo
qb.pulse_exp(exp = 'echo', n_avg = 6000, tmin = 0, tmax = 10e3, dt = 128, fit=True, check_mixers=False, flux=inst.get_ffl_bias()*1e3)

#%% T1

t_arr, I, Q,job,fitted_pars = qb.pulse_exp(exp='T1',n_avg=10000,tmin=0,tmax=50e3,dt=500,detuning=0)
#%% resonator ring down
qb.pulse_exp(exp = 'ringdown', n_avg = 30000, tmin = 4, tmax =5000, dt = 64, fit=True, check_mixers=False)

for i in np.round(np.linspace(91e6,100e6,10)):
    qb.update_value('qubit12_IF', i)
    qb.pulse_exp(exp = 'qubit_temp', n_avg = 2000, tmin = 40, tmax =240, dt = 4, fit=True, check_mixers=False, simulate=False, play_init_pi=True)
    



inst.turn_on_ffl_drive()

qb.pulse_exp(exp = 'cavity-reset', n_avg = 8000, tmin =0, tmax = 8e3, dt = 64, fit=True, check_mixers=False, simulate=True, ffl_len=0.6e3, with_ffl=True, amp_ffl_scale=0., amp_r_scale=0., detuning=0, flux=inst.get_ffl_bias()*1e3)

qb.pulse_exp(exp = 'cavity-cooling', n_avg = 6000, tmin =0, tmax = 10e3, dt = 128, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=1., amp_r_scale=0., detuning=0, flux=inst.get_ffl_bias()*1e3) 

qb.pulse_exp(exp = 'cavity-cooling-ramsey', n_avg = 15000, tmin =20, tmax = 8e3, dt = 128, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=1., amp_r_scale=0., detuning=-200e3, flux=inst.get_ffl_bias()*1e3) 

prog = qb.make_sequence(exp='ss', n_reps = 500, nIterations=1000)

datadict, job = qb.get_results(jobtype = prog,result_names=['i','I','Q','Iexc','Qexc'], showprogress=False, liveplot = False)

#plot, ax = pf.init_IQ_plot()
pf.plot_single_shot(datadict)

from qualang_tools.analysis.discriminator import two_state_discriminator

# Run an IQ blob experiment.

angle, threshold, fidelity, gg, ge, eg, ee = two_state_discriminator(datadict['I'], datadict['Q'], datadict['Iexc'], datadict['Qexc'], b_print=True, b_plot=True)

#%% diss T1
inst.turn_on_ffl_drive()

qb.update_value('rr_atten', 29)
qb.update_value('rr_pulse_len_in_clk', 1600)
t_arr, I, Q, job, fitted_pars = qb.pulse_exp(exp = 'dissT1', n_avg = 20e6, tmin = 16, tmax =200, dt = 4, fit=True, check_mixers=False)

#%% diss spec
I,Q,freqs,job = qb.qubit_spec(f_LO=3.2e9,amp_q_scaling=0.9,IF_min=50e3,IF_max=350e6,df=50e3,n_avg=1000,on_off=False,showprogress=True,savedata=False,check_mixers=False,)

#%% cavity ring down
qb.update_value('rr_atten', 29)
inst.set_attenuator(qb.pars['rr_atten'])
inst.set_ffl_LO(qb.pars['ffl_LO'])
inst.turn_on_ffl_drive()
I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=63e6,IF_max=93e6,df=0.1e6,n_avg=1000,savedata=True)
fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
qb.update_value('rr_freq', fc)
inst.set_attenuator(qb.pars['rr_atten'])
inst.get_attenuation()
dataDict, fig = measure_ringdown_drive_off(qb, tmax=4e3, dt=32,n_avg=100000)
inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
qb.update_value('ffl_freq', 2.15e9+50e6)
qb.update_value('ffl_LO', 2.1e9)
qb.update_value('ffl_IF', qb.pars['ffl_freq'] - qb.pars['ffl_LO'])
inst.set_ffl_LO(qb.pars['ffl_LO'])
qb.update_value('ffl_atten',12)
inst.set_ffl_attenuator(qb.pars['ffl_atten'])

inst.turn_on_ffl_drive()
I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=63e6,IF_max=93e6,df=0.1e6,n_avg=1000,savedata=True)
fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
qb.update_value('rr_freq', fc)
inst.set_attenuator(qb.pars['rr_atten'])
inst.get_attenuation()
qb.update_value('ffl_LO', 3.5e9)
inst.set_ffl_LO(qb.pars['ffl_LO'])
qb.update_value('ffl_IF', 50e6)
inst.set_ffl_LO(qb.pars['ffl_LO'])
dataDict, fig = measure_ringdown_drive_on(qb, amp_ffl_scale=0.9,tmax=1e3, dt=4, n_avg=100000)


dataDicton, fig, _ =measure_ringdown_drive_on(qb, amp_r_scale=1, tmin=0,  amp_ffl_scale=1.,tmax = 0.3e3,dt = 4,n_avg = 100000, flux=inst.get_ffl_bias()*1e3)
dataDictoff, fig= measure_ringdown_drive_off(qb, amp_r_scale=1 ,tmin=0, tmax = 2e3,dt = 32,n_avg = 10000,)

#%% qubit reset
dataDict, fig = measure_t1_ffl_off(qb, amp_r_scale=1,amp_ffl_scale=0.8, tmin = 0, tmax = 15e3, dt = 64, n_avg = 10000,)
# qb.update_value('ffl_freq', 3.8918e9)
#qb.update_value('ffl_freq', 3.9918e9)

inst.turn_on_ffl_drive()
dataDict, fig, pars = measure_t1_w_ffl(qb,amp_r_scale=1,amp_ffl_scale=1., tmin = 0, tmax = 100e3, dt = 512, n_avg = 2000)
dataDict, fig = measure_leakage_w_ffl(qb, amp_r_scale=1,amp_ffl_scale=1, tmin = 0, tmax =200e3, dt = 4096, n_avg = 1000,)
#for freq in [3.58e9]:


#%% T2 Values sweeping initial populating amp and wait time before first pi/2 pulse

t2vals=np.zeros((11,11))
for i, amp_r_scale in enumerate(np.linspace(0,0.35,11)):
    for j, amp_ffl_scale  in enumerate(np.linspace(0,0.36,11)):
        exp=qb.pulse_exp(exp = 'cavity-cooling', n_avg = 20000, tmin =0, tmax = 16e3, dt = 256, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=amp_ffl_scale, amp_r_scale=amp_r_scale, detuning=0, flux=inst.get_ffl_bias()*1e3) 
        #exp2=qb.pulse_exp(exp = 'cavity-reset', n_avg = 2000, tmin =0, tmax = 8e3, dt = 128, fit=True, check_mixers=False, simulate=False, ffl_len=ffl_len, with_ffl=True, amp_ffl_scale=amp_ffl_scale, amp_r_scale=0, detuning=500e3)
        t2vals[i,j]=exp[4][1]
       
import matplotlib.pyplot as plt       
fig, ax = plt.subplots(figsize=(6,6))
plot=ax.imshow(t2vals, interpolation='nearest', extent=[0,0.36,0,0.35], vmin=0.4, vmax=8, aspect="auto", origin="lower")
#ax.set_yticklabels([250, 500,1000,2000, 3000, 5000])
plt.xlabel("ffl amp")
plt.ylabel("rr amp")
plt.colorbar(plot)


datadict,fig=measure_ringdown_drive_on(qb, amp_ffl_scale=0.1, n_avg=2000, dt = 16,tmin=0, tmax=3e3,amp_r_scale=1, flux=inst.get_ffl_bias()*1e3)
datadictoff,fig=measure_ringdown_drive_off(qb, n_avg=2000, dt = 16,tmin=0, tmax=3e3,amp_r_scale=1,flux=inst.get_ffl_bias())

plt.plot(datadict["time"],datadict["Q"], label="drive on")
plt.plot(datadictoff["time"],datadictoff["Q"], label="drive off")
plt.xlabel("time")
plt.title("Q quad comparison")
plt.ylabel("Voltage")
plt.legend()
plt.show()

prog = qb.make_sequence(exp="opt_rr_freq", IFmin=202e6, IFmax=2045e5,df=100e3)

plt.plot(np.linspace(0,1,11),t2bare[5,:], label='without ffl')
plt.plot(np.linspace(0,1,11),t2vals[5,:], label= "with ffl")
plt.xlabel("rr amp scale")
plt.ylabel("T2 (in us)")
plt.title("Improvement in T2")
ax=plt.gca()
ax.set_ylim([0,4.6])
plt.legend()
plt.show()

T2woffl=np.zeros(13)
for i, scale in enumerate(np.linspace(0,0.13,13)):
    exp=qb.pulse_exp(exp = 'cavity-reset', n_avg = 3000, tmin = 100, tmax = 6e3, dt = 32, fit=True, check_mixers=False, simulate=True, ffl_len=200, with_ffl=True, amp_ffl_scale=0., amp_r_scale=scale, detuning=800e3)
    T2woffl[i]=exp[4][3]

T2wffl=np.zeros(13)
for i, scale in enumerate(np.linspace(0,0.13,13)):
    exp=qb.pulse_exp(exp = 'cavity-reset', n_avg = 3000, tmin = 0, tmax = 6e3, dt = 32, fit=True, check_mixers=False, simulate=True, ffl_len=200, with_ffl=True, amp_ffl_scale=0.13, amp_r_scale=scale, detuning=800e3)
    T2wffl[i]=exp[4][3]

for scale in np.linspace(0,1,11):
    qb.pulse_exp(exp = 'cavity-reset', n_avg = 2000, tmin =0, tmax = 60e3, dt = 512, fit=True, check_mixers=False, simulate=True, ffl_len=20, with_ffl=True, amp_ffl_scale=scale, amp_r_scale=1, detuning=0)
    
#%% 2 tone ffl sweep amp
from t1_ffl import measure_t1_w_ffl_qd_df
import matplotlib.pyplot as plt

mag_arr = []
amp_arr = np.linspace(0,1,51)

for i, amp in enumerate(amp_arr):
   
    mag = measure_t1_w_ffl_qd_df(qb, amp_ffl_scale=amp, amp_fflqc_scale=amp, n_avg = 2000, ffl_len = 1e3)
    mag_arr.append(mag)
    
plt.plot(amp_arr, mag_arr)
plt.xlabel("FFL amplitude")
plt.ylabel("Magnitude")
plt.plot()

#%% 2 tone ffl sweep flux
ring=np.zeros(30)
flux_arr = np.linspace(1.117,1.148,10) * 1e-3

for i, flux in enumerate(flux_arr):
    inst.set_ffl_bias(flux, step=10e-6, upper_bound = 10e-3, lower_bound = -10e-3)
    I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=30,IF_min=198e6,IF_max=203e6,df=50e3,n_avg=1000,savedata=True)
    fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
    qb.update_value('rr_freq', fc)
    qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'] )
    datadict,fig, fits=measure_ringdown_drive_on(qb, amp_ffl_scale=0.1, n_avg=10000, dt = 32,tmin=0, tmax=3e3,amp_r_scale=1, flux=inst.get_ffl_bias())
    ring[i]=fits[1]
    
    
plt.plot(flux_arr*1e3, mag_arr)
plt.xlabel("Flux (mA)")
plt.ylabel("Magnitude")
plt.plot()


measure_ringdown_drive_off(qb, n_avg=5000, dt = 16,tmin=0, tmax=2.5e3,amp_r_scale=1)


t2off=np.zeros(15)
t2on=np.zeros(15)
t2erroroff=np.zeros(15)
t2erroron=np.zeros(15)
for i, scale in enumerate(np.linspace(0.,0.35,15)):
    _,_,_,_,fitoff,erroroff=qb.pulse_exp(exp = 'cavity-cooling', n_avg = 15000, tmin =20, tmax = 8e3, dt = 64, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=0., amp_r_scale=scale, detuning=0, flux=inst.get_ffl_bias()*1e3) 
    _,_,_,_,fitson,erroron=qb.pulse_exp(exp = 'cavity-cooling', n_avg = 15000, tmin =20, tmax = 8e3, dt = 64, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=0.4, amp_r_scale=scale, detuning=0, flux=inst.get_ffl_bias()*1e3) 
    t2on[i]=fitson[1]
    t2erroron[i]=erroron[1]
    t2off[i]=fitsoff[1]
    t2erroroff[i]=erroroff[1]


qb.update_value('ffl_atten', 9)
inst.set_ffl_attenuator(qb.pars['ffl_atten'])
for flux in np.linspace(3e-3,6e-3,50):
    inst.set_ffl_bias(flux, step = 10e-6, lower_bound=-5e-3, upper_bound=6e-3)
    I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=30,IF_min=198e6,IF_max=210e6,df=50e3,n_avg=1000,savedata=True)
    fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
    qb.update_value('rr_freq', fc)
    qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'] )
    measure_ringdown_drive_on(qb, amp_ffl_scale=1, n_avg=2000, dt = 16,tmin=0, tmax=3e3,amp_r_scale=1, flux=flux*1e3)
    measure_t1_w_ffl(qb,amp_r_scale=1,amp_ffl_scale=1, tmin = 0, tmax = 60e3, dt = 1024, n_avg = 2000,flux=flux*1e3)
    
t1vals=np.zeros((10,51))
for i, flux in enumerate(np.linspace(0.7,1.2,15)*1e-3):
    inst.set_ffl_bias(flux, step = 10e-6, lower_bound=-5e-3, upper_bound=6e-3)
    I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=29,IF_min=196e6,IF_max=204e6,df=50e3,n_avg=1500,savedata=True)
    fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
    qb.update_value('rr_freq', fc)
    qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'] )
    I,Q,freqs,job = qb.run_scan(df = 2e6, n_avg =6000, element='fflqb', chunksize = 350e6,  lo_min = 2e9, lo_max = 10e9, amp_q_scaling = 1, saturation_dur = 25e3, showprogress=True, res_ringdown_time = int(40e3), check_mixers=True, flux=inst.get_ffl_bias()*1e3)
    #qb.pulse_exp(exp='T1', check_mixers = False, n_avg=2000,dt=1024,tmin=0,tmax=60e3,detuning=0e3, amp_ffl_scale=0, flux=inst.get_ffl_bias()*1e3)
    # for j, ffl_IF  in enumerate(np.linspace(0,300,51)):
    #     qb.update_value('ffl_IF', ffl_IF)
    #     exp=measure_t1_w_ffl(qb,amp_r_scale=1,amp_ffl_scale=0.75, tmin = 0, tmax = 60e3, dt = 1024, n_avg = 2000,)
    #     if exp[2][0]>0:
    #         t1vals[i,j]=exp[2][1]
    #     else:
    #         t1vals[i,j]=0
            
fig, ax = plt.subplots(figsize=(6,6))
plot=ax.imshow(t1vals, interpolation='nearest', extent=[0,300,-200,-1000], vmin=10, vmax=25, aspect="auto",)
ax.set_yticklabels(np.linspace(-0.200,-1,10))
plt.xlabel("4.95 GHz+ IF (MHz)")
plt.ylabel("flux (mA)")
plt.colorbar(plot)

from sequence import *
host='10.71.0.56'
port='9510'
prog = qb.make_sequence(exp='ss_cav', n_reps = 1000, nIterations=30)
qmm = QuantumMachinesManager(host=host, port=port)
job = qmm.simulate(config=qb.config, program=prog, simulate=SimulationConfig(duration=3000))
job.get_simulated_samples().con1.plot()

#%paper plots%

dataDicton, fig=measure_ringdown_drive_on(qb, amp_r_scale=1, tmin=0,  amp_ffl_scale=0.15,tmax = 4e3,dt = 16,n_avg = 2000,)
dataDictoff, fig=measure_ringdown_drive_off(qb, amp_r_scale=1 ,tmin=96, tmax = 4e3,dt = 16,n_avg = 100000,)
inst.turn_on_ffl_drive()
qb.pulse_exp(exp = 'cavity-cooling', n_avg = 40000, tmin =20, tmax = 10e3, dt = 64, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=0., amp_r_scale=0., detuning=0, flux=inst.get_ffl_bias()*1e3) 
qb.pulse_exp(exp = 'cavity-cooling', n_avg = 25000, tmin =20, tmax = 5e3, dt = 32, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=0., amp_r_scale=0.25, detuning=0, flux=inst.get_ffl_bias()*1e3) 
qb.pulse_exp(exp = 'cavity-cooling', n_avg = 30000, tmin =20, tmax = 10e3, dt = 64, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=0.35, amp_r_scale=0.25, detuning=0, flux=inst.get_ffl_bias()*1e3) 
qb.pulse_exp(exp = 'cavity-cooling', n_avg = 10000, tmin =20, tmax = 16e3, dt = 128, fit=True, check_mixers=False, simulate=True, with_ffl=True, amp_ffl_scale=0.35, amp_r_scale=0.35, detuning=0, flux=inst.get_ffl_bias()*1e3) 

inst.turn_off_ffl_drive()
qb.pulse_exp(exp = 'cavity-cooling-ramsey', n_avg = 30000, tmin =20, tmax = 10e3, dt = 128, fit=True, check_mixers=False, simulate=False, with_ffl=True, amp_ffl_scale=0., amp_r_scale=0., detuning=-200e3, flux=inst.get_ffl_bias()*1e3) 
inst.turn_on_ffl_drive()

det=np.zeros(21)
coh=np.zeros(21)
for i, amp in enumerate(np.linspace(0.2,0.03,21)):
    fitted_pars = qb.pulse_exp(exp='ramsey_chi', check_mixers = False, n_avg=6000,dt=256,tmin=0,tmax=12e3,detuning=0.1e6,amp_r_scale=amp,simulate=False)
    det[i]=fitted_pars[1]
    coh[i]=fitted_pars[3]
qb.pulse_exp(exp = 'cavity-cooling-ramsey', n_avg = 40000, tmin =20, tmax = 10e3, dt = 128, fit=True, check_mixers=False, simulate=False, with_ffl=True, amp_ffl_scale=0., amp_r_scale=0., detuning=-200e3, flux=inst.get_ffl_bias()*1e3)     
    
    
for i in np.linspace(200.20e6, 200.40e6, 11):
    qb.update_value('rr_IF', 200.36e6)
    prog = qb.make_sequence(exp='ss', n_reps = 300,nIterations=2000)
    datadict, job = qb.get_results(jobtype = prog,result_names=['i', 'I','Q','Iexc','Qexc'], showprogress=False, liveplot = False)
    angle, threshold, fidelity, gg, ge, eg, ee = two_state_discriminator(datadict['I'], datadict['Q'], datadict['Iexc'], datadict['Qexc'], b_print=True, b_plot=True)
    #pf.plot_single_shot(datadict)
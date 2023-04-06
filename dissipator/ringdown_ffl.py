# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 18:29:56 2023

@author: lfl
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 18:14:15 2023

@author: lfl
"""
import timeit
from qubit import *
import instrument_init as inst
import plot_functions as pf
import h5py
from datetime import datetime
import os

qb = qubit('logical')
qb.update_value('ffl_freq', 2.85e9)
qb.update_value('ffl_LO', 2.8e9)
qb.update_value('ffl_IF', 0.05e9)
qb.update_value('amp_ffl', 0.45)
inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
inst.set_rr_LO(qb.pars['rr_LO']) # turn on
bOptimizeMixer = False

from instrument_init import init_sa, init_sa_by_serial_number

if bOptimizeMixer:
    # mixer optimization
    ref_H = 0
    ref_L = -30
    qb.play_pulses()
    sa = init_sa_by_serial_number(20234492)
    qb_lo_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO'],reference=ref_H,config=True,plot=True)
    qb_im_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO']-qb.pars['ffl_IF'],reference=ref_H,config=True,plot=True)
    qb_on_power = qb.get_power(sa, freq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
    
    '''Optimize FFL mixer'''
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
    # qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
    # qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element = 'ffl')
    
    sa_close_device(sa)
    sa = init_sa_by_serial_number(20234229)
    set_attenuator(0)
    get_attenuation()
    rr_lo_leakage = qb.get_power(sa, freq=qb.pars['rr_LO'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
    rr_im_leakage = qb.get_power(sa, freq=qb.pars['rr_LO']-qb.pars['rr_IF'],span = 1e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
    rr_on_power = qb.get_power(sa, freq=qb.pars['rr_LO']+qb.pars['rr_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

    # do a coarse sweep to minimize LO leakage
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element='rr')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine',element='rr')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element='rr')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='intermediate',reference = ref_H, element='rr')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element='rr')
    sa_close_device(sa)
    
    qb.write_pars()

def measure_ringdown_drive_on(amp_r_scale=1, 
                              amp_ffl_scale=1,
                              tmin = 0,
                              tmax = 2e3,
                              dt = 16,
                              n_avg = 1000,
                              ):
    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    
    inst.set_attenuator(attenuation=qb.pars['rr_atten'])
    inst.turn_on_ffl_drive()
    resettime_clk= clk(qb.pars['rr_resettime'])
    with program() as prog:
        update_frequency('rr', (qb.pars['rr_freq']-qb.pars['rr_LO'])) 
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO'])) 
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    play("readout"*amp(amp_r_scale), "rr")
                    measure("void", "rr", None,*qb.res_demod(I,Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("readout"*amp(amp_r_scale), "rr")
                    align("ffl", "rr")
                    play('const'*amp(amp_ffl_scale), "ffl", duration=t)
                    align("ffl", "rr")
                    measure("void", "rr", None,*qb.res_demod(I,Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
    
        with stream_processing():
            I_stream.buffer(len(t_arr)).average().save("I")
            Q_stream.buffer(len(t_arr)).average().save("Q")
            n_stream.save('n')

    datadict, job = qb.get_results(jobtype = prog, result_names = ["I", "Q"], n_total=n_avg, notify = False)
    t_arr = np.array(t_arr)*4/1e3
    I = np.array(datadict["I"])
    Q = np.array(datadict["Q"])
    ydata = np.abs(I+1j*Q)

    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='ringdown',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='ringdown_drive=ON',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                 qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                             'amp_ffl_scale': amp_ffl_scale,
                             'tmin': 16,
                             'tmax': 2e3,
                             'dt': 16,
                             'n_avg': 1000,},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }
    return dataDict, fig

def measure_ringdown_drive_off(amp_r_scale=1,
                               tmin =  0,
                               tmax= 2e3,
                               dt= 16,
                               n_avg= 1000,
                               ):
    inst.turn_off_ffl_drive()
    inst.set_attenuator(attenuation=qb.pars['rr_atten'])
    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    resettime_clk= clk(qb.pars['rr_resettime'])
    
    with program() as prog:
        update_frequency('rr', (qb.pars['rr_freq']-qb.pars['rr_LO'])) 
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO']))
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])

        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)

        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    play("readout" * amp(amp_r_scale), "rr")
                    measure("void", "rr", None,*qb.res_demod(I,Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("readout" * amp(amp_r_scale), "rr")
                    wait(t, 'rr')
                    measure("void", "rr", None,*qb.res_demod(I,Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)

        with stream_processing():
            I_stream.buffer(len(t_arr)).average().save("I")
            Q_stream.buffer(len(t_arr)).average().save("Q")
            n_stream.save('n')

    datadict, job = qb.get_results(jobtype = prog, result_names = ["I", "Q"], n_total=n_avg, notify = False)
    t_arr = np.array(t_arr)*4/1e3
    I = np.array(datadict["I"])
    Q = np.array(datadict["Q"])
    ydata = np.abs(I+1j*Q)
    
    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='ringdown_off',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='ringdown_drive=off',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                 qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
    
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                             'ffl_drive_output': False ,
                             'tmin': 16,
                             'tmax': 2e3,
                             'dt': 16,
                             'n_avg': 1000,},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }

    
    return dataDict, fig

def save_datadict_to_fgroup(f, name, datadict):
    subgroup = f.create_group(name)
    dset_i = subgroup.create_dataset('I', data=datadict['I'])
    dset_q = subgroup.create_dataset('Q', data=datadict['Q'])
    dset_t = subgroup.create_dataset('t', data=datadict['time'])
    for key in datadict['metadata'].keys():
        subgroup.attrs[key] = datadict['metadata'][key]
    print(f'write dataset to {name}')
    
    


#%% doing sweeps
def sweep_powers(ffl_freq = 6.6e9,rr_atten = 23, n_avg=4000):
    qb.update_value('ffl_freq', ffl_freq)
    ffl_scales = np.round(np.linspace(0.0,1.0,20),2)
    # rr_scales = np.round(np.linspace(0.1,1,10), 2)
    rr_scales = [1]
    avgs = [n_avg/r for r in rr_scales]
    rr_atten = qb.pars['rr_atten']
    rrLen = qb.pars['rr_pulse_len_in_clk']
    inst.set_attenuator(rr_atten)
    
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    filename = f'ringdown_sweepPowers_fflFreq={(ffl_freq)/1e9:.2f}GHz_DA={rr_atten}dB_rrLen={rrLen}clks_navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_on = hf.create_group(f'sweep_ffl_amp_{timestamp}')
        for i,sr in enumerate(rr_scales):
            dataDict, fig = measure_ringdown_drive_off(amp_r_scale=sr, tmax=1e3, dt=4, n_avg=n_avg)
            save_datadict_to_fgroup(g_on, f'ro_scale = {sr}, ffl_off', dataDict)
          
            for j,sf in enumerate(ffl_scales):
                dataDict, fig = measure_ringdown_drive_on(amp_r_scale=sr, amp_ffl_scale=sf, n_avg=n_avg, dt = 4, tmax=1e3)
                save_datadict_to_fgroup(g_on, f'ro_scale = {sr}, ffl_scale = {sf}', dataDict)
                
def main():
    # save data
    qb.update_value('rr_pulse_len_in_clk',int(500)) # default 500
    qb.update_value('rr_atten', 35) # default = 23
    device = 'diss08_07A'
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\ringdown'
    
    start = timeit.default_timer()   
    # LO_freqs = np.linspace(5.5e9, 7.5e9, 5)
    for rr_atten in [35]:
        # do res spec
        inst.turn_off_ffl_drive()
        I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=rr_atten,IF_min=30e6,IF_max=60e6,df=0.1e6,n_avg=1000,savedata=True)
        fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
        qb.update_value('rr_freq', fc)
        dataDict = measure_ringdown_drive_off(amp_r_scale=1, tmax=4e3, dt=16, n_avg=10000)
        inst.turn_on_ffl_drive()
        sweep_powers(ffl_freq = qb.pars['ffl_freq'],rr_atten = rr_atten,n_avg=100000)
        stop = timeit.default_timer()
        print('Time: ', stop - start)  

    
if __name__ == "__main__":
    main()
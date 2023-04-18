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
from instrument_init import init_sa, init_sa_by_serial_number

def optimize_mixer(sa, qb, element='rr', cal='LO'):
    ref_H = 0
    ref_L = -45
    qb.play_pulses(element)
    if element == 'rr':
        inst.set_attenuator(0)
        inst.get_attenuation()
    qb_lo_leakage = qb.get_power(sa, freq=qb.pars[f'{element}_LO'],reference=ref_H,config=True,plot=True)
    qb_im_leakage = qb.get_power(sa, freq=qb.pars[f'{element}_LO']-qb.pars[f'{element}_IF'],reference=ref_H,config=True,plot=True)
    qb_on_power = qb.get_power(sa, freq=qb.pars[f'{element}_LO']+qb.pars[f'{element}_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
    
    '''Optimize FFL mixer'''
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='coarse', element = element)
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='intermediate', element = element)
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='intermediate', element = element)
    # qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_L, mode='fine', element = element)
    
    qb.write_pars()

def measure_ringdown_drive_on(qb,
                              amp_r_scale=1, 
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
    fig = pf.plot_data(t_arr,ydata,sequence='ringdown_drive=ON',fitted_pars=fitted_pars,nAverages=n_avg, 
                 qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1)
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

def measure_ringdown_drive_off(qb,
                               amp_r_scale=1,
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
    fig = pf.plot_data(t_arr,ydata,sequence='ringdown_drive=off',fitted_pars=fitted_pars,nAverages=n_avg, 
                 qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1)
    
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
    
#%% ffl punchout
def ffl_punchout():
    qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['rr_freq'])
    qb.update_value('ffl_LO', 3.87e9)
    qb.update_value('ffl_IF', qb.pars['ffl_freq'] - qb.pars['ffl_LO'])
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    stepsize = 2
    atten_list = np.arange(30,2,-2)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    data = []
    for atten in atten_list:
        inst.set_ffl_attenuator(atten)
        I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=63e6,IF_max=93e6,df=0.1e6,n_avg=1000,savedata=True)
        data.append(np.abs(I + 1j * Q))
    im = ax.imshow(data, aspect='auto',origin='lower',extent=(freqs[0]/1e9, freqs[-1]/1e9,atten_list[0], atten_list[-1]),
            interpolation=None, cmap='RdBu')
    ax.legend()
    plt.show()
    return fig



#%% doing sweeps
def sweep_powers(qb,ffl_freq = 6.6e9,rr_atten = 23, n_avg=4000):
    rr_atten = qb.pars['rr_atten']
    rrLen = qb.pars['rr_pulse_len_in_clk']
    ffl_scales = [round(0.0 + 0.1* n,1) for n in range(10)]
   
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    filename = f'ringdown_sweepPowers_flux=70uA_fflFreq={(ffl_freq)/1e9:.2f}GHz_DA={rr_atten}dB_fDA={qb.pars["ffl_atten"]}dB_rrLen={rrLen}clks_navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_on = hf.create_group(f'sweep_ffl_amp_{timestamp}')
        
        dataDict, fig = measure_ringdown_drive_off(qb,tmax=2e3, dt=16, n_avg=n_avg)
        save_datadict_to_fgroup(g_on, f' ffl_off', dataDict)
      
        for j,sf in enumerate(ffl_scales):
            dataDict, fig = measure_ringdown_drive_on(qb,amp_ffl_scale=sf, n_avg=n_avg, dt = 16, tmax=2e3)
            save_datadict_to_fgroup(g_on, f'ffl_scale = {sf}', dataDict)
                
def main():
    qb = qubit('diss09')
    qb.update_value('diss_freq', 9.5e9)
    qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['rr_freq'])
    
    qb.update_value('ffl_IF', 250e6)
    qb.update_value('ffl_LO', qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
    qb.update_value('amp_ffl', 0.45)
    qb.update_value('ffl_atten', 30)
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    inst.set_rr_LO(qb.pars['rr_LO']) # turn on
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    # sweep parameters
    IF_min = 50e6
    IF_max = 250e6
    stepsize = 50e6
    # IF_list = [int(qb.pars['ffl_IF'] + stepsize * n) for n in range(-4,3)]
    IF_list = [qb.pars['ffl_IF']]
    bOptimizeFFLMixer = True
    bOptimizeRRMixer = True
    # save data
    device = 'diss09'
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\ringdown'
    
    start = timeit.default_timer() 
    # readout mixer optimization
    sa = inst.init_sa()
    if bOptimizeRRMixer:
        qb.play_pulses(element='rr')
        inst.set_attenuator(0)
        optimize_mixer(sa, qb, element='rr',cal='LO')
        optimize_mixer(sa, qb, element='rr',cal='SB')
        # sa_close_device(sa)

    # do res spec
    inst.turn_off_ffl_drive()
    I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=63e6,IF_max=93e6,df=0.1e6,n_avg=1000,savedata=True)
    fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
    qb.update_value('rr_freq', fc)
    
    # change the following line to sweep LO frequencies
    lo_list = [qb.pars['ffl_LO'] - 350e6 * n for n in [0]]
    for ffl_lo in lo_list:
        inst.set_ffl_LO(ffl_lo)

        if bOptimizeFFLMixer:
            # sa = inst.init_sa()
            qb.play_pulses(element='ffl')
            optimize_mixer(sa, qb, element='ffl',cal='LO')
            optimize_mixer(sa, qb, element='ffl',cal='SB')
            # sa_close_device(sa)
       
        I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=63e6,IF_max=93e6,df=0.1e6,n_avg=1000,savedata=True)
        fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
        qb.update_value('rr_freq', fc)
        for ffl_IF in IF_list:
            qb.update_value('ffl_IF', ffl_IF)
            qb.update_value('ffl_freq', ffl_lo + ffl_IF)
            if bOptimizeFFLMixer:
                # ffl mixer calibration
                qb.play_pulses('ffl')
                optimize_mixer(sa, qb, element='ffl',cal='SB')
                # sa_close_device(sa)
            sweep_powers(qb,ffl_freq = qb.pars['ffl_freq'],rr_atten = qb.pars['rr_atten'],n_avg=100000)
    stop = timeit.default_timer()
    sa_close_device(sa)
    print('Time: ', stop - start)  

    
if __name__ == "__main__":
    main()
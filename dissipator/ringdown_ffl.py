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
from dissipator import *
from sequence import *
import instrument_init as inst
import plot_functions as pf
import h5py
from datetime import datetime
import os
from instrument_init import init_sa, init_sa_by_serial_number
from fflCalibration import acquire_rr_spec_background
from Utilities import clk
from g_e_f_T1 import measure_leakage_w_ffl
def measure_ringdown_drive_on(qb,
                              amp_r_scale=1, 
                              amp_ffl_scale=1,
                              tmin = 0,
                              tmax = 2e3,
                              dt = 16,
                              n_avg = 1000,
                              flux=0):

    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    
    inst.set_attenuator(attenuation=qb.pars['rr_atten'])
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    inst.turn_on_ffl_drive()
    
    seq = sequence('ringdown_drive_on',n_avg=n_avg, amp_r_scale=amp_r_scale, amp_ffl_scale=amp_ffl_scale)
    prog = seq.make_sequence(qb,tmin=tmin, tmax=tmax, dt=dt)

    datadict, job = qb.get_results(jobtype = prog, result_names = ["I", "Q"], n_total=n_avg, notify = False)
    t_arr = np.array(t_arr)/1e3
    I = np.array(datadict["I"])
    Q = np.array(datadict["Q"])
    ydata = np.abs(I+1j*Q)
    #ydata=-Q
    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='ringdown',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='ringdown_drive=ON',fitted_pars=fitted_pars,nAverages=n_avg, 
                 qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1, flux=flux, amp_ffl_scale=amp_ffl_scale, rr_atten=qb.pars['rr_atten'], ffl_atten=qb.pars['ffl_atten'], error=error)
    dataDict = {'metadata': {'ffl_freq': qb.pars['ffl_freq'],
                             'ffl_atten': qb.pars['ffl_atten'],
                             'ffl_LO': qb.pars['ffl_LO'],
                             'ffl_IF': qb.pars['ffl_IF'],
                             'rr_freq': qb.pars['rr_freq'],
                             'rr_LO': qb.pars['rr_LO'],
                             'rr_IF': qb.pars['rr_IF'],
                             'amp_ffl': qb.pars['amp_ffl'],
                             'amp_r': qb.pars['amp_r'],
                             'amp_r_scale': amp_r_scale,
                             'amp_ffl_scale': amp_ffl_scale,
                             'tmin': 16,
                             'tmax': 2e3,
                             'dt': 16,
                             'n_avg': 1000,
                             'report': str(job.execution_report()),
                             },
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }
    rrLen = qb.pars['rr_pulse_len_in_clk']
    expName= 'cavity_ringdown'
    filename = f'{expName}_flux_ffl_len={flux*1e3}mA_fflFreq={(qb.pars["ffl_freq"])/1e9:.2f}GHz_DA={qb.pars["rr_atten"]}dB_fDA={qb.pars["ffl_atten"]}dB_rrLen={rrLen}clks_navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)
    #inst.set_attenuator(qb.pars['rr_atten'])
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        flux_key = f'flux={flux*1e6:.1f}uA'
        freq_key = f'ffl_freq={qb.pars["ffl_freq"]/1e9:.5f}GHz'
        #ffl_len_key= f'ffl_len_key={ffl_len:.1f}us'
        data_key = '_'.join([flux_key,freq_key,timestamp])
        g_on = hf.create_group(data_key)
        save_datadict_to_fgroup(g_on, f'ffl amp = {amp_ffl_scale:.3f}', dataDict)
    return dataDict, fig, fitted_pars

def measure_ringdown_drive_off(qb,
                               amp_r_scale=1,
                               tmin =  0,
                               tmax= 2e3,
                               dt= 16,
                               n_avg= 1000,flux=0
                               ):
    inst.turn_off_ffl_drive()
    inst.set_attenuator(attenuation=qb.pars['rr_atten'])
    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    amp_ffl_scale=0
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    resettime_clk= clk(qb.pars['rr_resettime'])
    
    with program() as prog:
        update_frequency('rr', (qb.pars['rr_freq']-qb.pars['rr_LO'])) 
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO']))
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])

        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)

        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_each_(t, t_arr):
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
    ydata = (np.abs(I+1j*Q))**2
    #ydata=I**2
    
    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='ringdown_off',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='ringdown_drive=off',fitted_pars=fitted_pars,nAverages=n_avg, 
                 qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1)
    
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                             'ffl_drive_output': False ,
                             'tmin': 16,
                             'tmax': 2e3,
                             'dt': 16,
                             'n_avg': 1000,
                             'report': str(job.execution_report())},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }

    rrLen = qb.pars['rr_pulse_len_in_clk']
    expName= 'cavity_ringdown'
    filename = f'{expName}_flux_ffl_len={flux*1e3}mA_fflFreq={(qb.pars["ffl_freq"])/1e9:.2f}GHz_DA={qb.pars["rr_atten"]}dB_fDA={qb.pars["ffl_atten"]}dB_rrLen={rrLen}clks_navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        flux_key = f'flux={flux*1e6:.1f}uA'
        freq_key = f'ffl_freq={qb.pars["ffl_freq"]/1e9:.5f}GHz'
        #ffl_len_key= f'ffl_len_key={ffl_len:.1f}us'
        data_key = '_'.join([flux_key,freq_key,timestamp])
        g_on = hf.create_group(data_key)
        save_datadict_to_fgroup(g_on, f'ffl amp = {amp_ffl_scale:.3f}', dataDict)
    return dataDict, fig

def save_datadict_to_fgroup(f, name, datadict):
    subgroup = f.create_group(name)
    for k in datadict.keys():
        if k != 'metadata':
            subgroup.create_dataset(k, data=datadict[k])
   
    for key in datadict['metadata'].keys():
        subgroup.attrs[key] = datadict['metadata'][key]
    print(f'write dataset to {name}')

    
#%% ffl punchout
def ffl_punchout(qb,stepsize = 2, save=True):
    if save:
        today = datetime.today()
        sDate =  today.strftime("%Y%m%d")
        saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\punchout'
    
        if not os.path.exists(saveDir):
            Path(saveDir).mkdir(parents=True, exist_ok=True)
        filename = f'fflPunchout_fflFreq={(qb.pars["ffl_freq"])/1e9:.2f}GHz_DA={rr_atten}dB_amp_r_scale={amp_r_scale}'
        index = get_index_for_filename(saveDir, filename)
    qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['rr_freq'])
    qb.update_value('ffl_IF', 350e6)
    qb.update_value('ffl_LO', qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
    
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    atten_list = np.arange(30,-2,-stepsize)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    data = []
    for atten in atten_list:
        inst.set_ffl_attenuator(atten)
        I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=35e6,IF_max=65e6,df=0.1e6,n_avg=1000,savedata=True)
        data.append(np.abs(I + 1j * Q))
    im = ax.imshow(data, aspect='auto',origin='lower',extent=(freqs[0]/1e9, freqs[-1]/1e9,atten_list[0], atten_list[-1]),
            interpolation=None, cmap='RdBu')
    ax.legend()
    plt.show()
    return fig



#%% doing sweeps
def sweep_powers(qb, 
                 sa = None,
                 element='ffl', 
                 n_avg=4000, 
                 check_mixers=False,
                 amp_r_scale=1, 
                 amp_ffl_scale=1,
                 bcheckDriveOff = False,
                 **kwargs):
    rr_atten = qb.pars['rr_atten']
    rrLen = qb.pars['rr_pulse_len_in_clk']
    if element == 'rr':
        var_scales = [round(0.0 + 0.1* n,1) for n in range(1,11)]
        var_scales.reverse()
        avg_list = [n_avg * n for n in range(1,11) ]
    elif element == 'ffl':
        var_scales = [9,11,13,15,17,20,23,26,30,34,38]
        #var_scales=[0]
    if check_mixers:
        if sa is None:
            sa = inst.init_sa()
        # qb.play_pulses(element='ffl')
        # qb.optimize_mixer(sa,element='ffl',cal='LO')
        # qb.optimize_mixer(sa, element='ffl',cal='SB')
        mixer_data_on = qb.mixer_powers(sa, element='ffl')
        mixer_data_off = qb.mixer_powers(sa, element='ffl', switch='off')
        # sa_close_device(sa)
     
    #flux = 0
    min_atten=qb.pars['ffl_atten']
    flux=inst.get_ffl_bias()
    flux_key = f'flux={flux*1e6:.1f}uA'
    freq_key = f'ffl_freq={qb.pars["ffl_freq"]/1e9:.5f}GHz'
    if 'filename' in kwargs.keys():
        saveDir = kwargs.get('saveDir')
        filename = kwargs.get('filename')
        index = kwargs.get('index')
        timestamp = kwargs.get('timestamp')
        mode = 'a'
    else:
        # save data
        today = datetime.today()
        sDate = today.strftime("%Y%m%d")
        saveDir = f'G:\\Shared drives\\CavityCooling\data\\{qb.device_name}\\{sDate}\\ringdown'
        
        if not os.path.exists(saveDir):
            Path(saveDir).mkdir(parents=True, exist_ok=True)
      
        filename = f'ringdown_sweep{element}Powers_flux={flux*1e6}uA_fflFreq={(qb.pars["ffl_freq"])/1e9:.2f}GHz_DA={rr_atten}dB_fDA={min_atten}_rrLen={rrLen}clks_amp_r_scale={amp_r_scale}=navg={n_avg}_checkMixer={str(check_mixers)}'
        index = get_index_for_filename(saveDir, filename)
        mode = 'w'
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5',mode) as hf:
        g_data = hf[f'ringdown_{timestamp}']
        if 'data_key' in kwargs.keys():
            data_key = kwargs.get('data_key')
        else:
            now = datetime.now()
            timestamp = now.strftime("%H:%M:%S")
            data_key = '_'.join([flux_key,freq_key,timestamp])
            if flux_key not in data_key:
                raise ValueError(f'flux mismatch: set value {flux_key} does not match data key {data_key}')
            if freq_key not in data_key:
                raise ValueError(f'ffl_freq mismatch: set value {freq_key} does not match data key {data_key}')
       
        g_subdata = g_data.create_group(data_key)
        
        if bcheckDriveOff:
            
            dataDict, fig = measure_ringdown_drive_off(qb,tmax=2.5e3, dt=16, n_avg=n_avg, amp_r_scale = amp_r_scale)
            save_datadict_to_fgroup(g_subdata, f' ffl_off', dataDict)
      
        for j,value in enumerate(var_scales):
            if element == 'rr':
                #dataDict, fig = measure_ringdown_drive_on(qb,amp_ffl_scale=amp_ffl_scale, n_avg=avg_list[j], dt = 8, tmin=96,tmax=0.8e3,amp_r_scale=value)
                pass
            elif element =='ffl':
                qb.update_value('ffl_atten', value)
                inst.set_ffl_attenuator(qb.pars['ffl_atten'])
                dataDict, fig, _ = measure_ringdown_drive_on(qb,amp_ffl_scale=1., n_avg=n_avg, dt = 8, tmin=0, tmax=1.8e3,amp_r_scale=amp_r_scale, flux=flux*1e3)
                #inst.turn_on_ffl_drive()
                #qb.pulse_exp(exp = 'cavity-cooling', n_avg = 6000, tmin =0, tmax = 10e3, dt = 128, fit=True, check_mixers=False, simulate=False, with_ffl=True, amp_ffl_scale=1., amp_r_scale=0., detuning=0, flux=inst.get_ffl_bias()*1e3) 
                save_datadict_to_fgroup(g_subdata, f'{element}_atten = {value}', dataDict)
        #qb.update_value('ffl_atten', min_atten)
        
        # meta data for mixer calibration

def sweep_frequency(qb, 
                 sa = None,
                 element='ffl', 
                 stepsize = 50e6,
                 fmin = 3.4e9,
                 fmax = 4.6e9,
                 n_avg=4000, 
                 check_mixers=False,
                 amp_r_scale=1, 
                 amp_ffl_scale=1,
                 bcheckDriveOff = False):
    rr_atten = qb.pars['rr_atten']
    rrLen = qb.pars['rr_pulse_len_in_clk']
    # qb.update_value('amp_ffl', 0.4)
    qb.update_value('ffl_atten', 10)
    qb.update_value('ffl_IF', 350e6)
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    inst.set_rr_LO(qb.pars['rr_LO']) # turn on
    # readout spec
    #I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=35e6,IF_max=65e6,df=0.1e6,n_avg=1000,savedata=True)
    # fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
    # qb.update_value('rr_freq', fc)
    # save data
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{qb.device_name}\\{sDate}\\ringdown'
    
    num_of_points = int((fmax - fmin)/stepsize) + 1
    freqs_list = [fmin + n * stepsize for n in range(num_of_points)]
    
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    flux = inst.get_flux_bias()
    #flux=0
    filename = f'ringdown_sweep{element}Freq_flux={int(flux*1e6)}uA_ffl_IF={qb.pars["ffl_IF"]/1e6}MHz_amp_ffl_scale={amp_ffl_scale}_DA={rr_atten}dB_fDA={qb.pars["ffl_atten"]}dB_rrLen={rrLen}clks_amp_r_scale={amp_r_scale}=navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)

    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_on = hf.create_group(f'sweep_{element}_freq_{timestamp}')
        
        if bcheckDriveOff:
            dataDict, fig = measure_ringdown_drive_off(qb,tmax=2e3, dt=16, n_avg=n_avg, amp_r_scale = amp_r_scale)
            save_datadict_to_fgroup(g_on, f'ffl_off', dataDict)
      
        for j,freq in enumerate(freqs_list):
            qb.update_value('ffl_freq', freq)
            qb.update_value('ffl_LO', qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
            inst.set_ffl_LO(qb.pars['ffl_LO'])
            if check_mixers:
                sa = inst.init_sa()
                qb.optimize_mixer(sa,element='ffl',cal='LO')
                # qb.optimize_mixer(sa, element='ffl',cal='SB')
                mixer_data_on = qb.mixer_powers(sa, element='ffl')
                sa_close_device(sa)
            print(f'ffl_IF={qb.pars["ffl_IF"]} Hz')
            dataDict, fig = measure_ringdown_drive_on(qb,amp_ffl_scale=amp_ffl_scale, n_avg=n_avg, dt = 16, tmax=2e3,amp_r_scale=amp_r_scale)
            # add metadata
            dataDict['metadata']['ffl_mixer_on'] = mixer_data_on 
            save_datadict_to_fgroup(g_on, f'{element}_freq = {freq}', dataDict)
        
        # meta data for mixer calibration
    
        
def sweep_freq_power(qb, sa = None,  
                     stepsize = 50e6,
                     fmin = 1e9,
                     fmax = 6e9,
                     n_avg=1000, 
                     bOptimizeFFLMixer=True,test=False, **kwargs):
   
    num_of_points = int((fmax - fmin)/stepsize) + 1
    freqs_list = [fmin + n * stepsize for n in range(num_of_points)]
    if 'filename' in kwargs:
        saveDir = kwargs.get('saveDir')
        filename = kwargs.get('filename')
        index = kwargs.get('index')
        timestamp = kwargs.get('timestamp')
        mode = 'a'
    else:
        # save data
        today = datetime.today()
        sDate =  today.strftime("%Y%m%d")
        saveDir = f'G:\\Shared drives\\CavityCooling\data\\{qb.device_name}\\{sDate}\\ringdown'
        
        if not os.path.exists(saveDir):
            Path(saveDir).mkdir(parents=True, exist_ok=True)
      
        filename = f'ringdown_sweep{element}Freq_flux={int(flux*1e6)}uA_ffl_IF={qb.pars["ffl_IF"]/1e6}MHz_amp_ffl_scale={amp_ffl_scale}_DA={rr_atten}dB_fDA={qb.pars["ffl_atten"]}dB_rrLen={rrLen}clks_amp_r_scale={amp_r_scale}=navg={n_avg}'
        index = get_index_for_filename(saveDir, filename)
        mode = 'w'
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5',mode) as hf:
        g_data = hf[f'ringdown_{timestamp}']
        g_calib = hf[f'calib_{timestamp}']
        flux = inst.get_ffl_bias()
        #flux=0
        flux_key = f'flux={flux*1e6:.1f}uA'
        #qb.update_value('ffl_mixer_offsets', [0.0,0.0])
        #qb.update_value('ffl_mixer_imbalance', [0.0,0.0])
        for freq in tqdm(freqs_list):
            qb.update_value('ffl_freq', freq)
            qb.update_value('ffl_LO', qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
            inst.set_ffl_attenuator(qb.pars['ffl_atten'])
            inst.set_rr_LO(qb.pars['rr_LO']) # turn on
            inst.set_ffl_LO(qb.pars['ffl_LO'])
            
            freq_key = f'ffl_freq={qb.pars["ffl_freq"]/1e9:.5f}GHz'
            now = datetime.now()
            timestamp = now.strftime("%H:%M:%S")
            data_key = '_'.join([flux_key,freq_key,timestamp])
            g_subcalib = g_calib.create_group(data_key)
            # do res spec
            # inst.turn_off_ffl_drive()
            # I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=141e6,IF_max=148e6,df=50e3,n_avg=5000,savedata=True, fit=True)
            # fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
            # qb.update_value('rr_freq', fc)
            # qb.update_value('rr_IF', qb.pars['rr_freq']-qb.pars['rr_LO'])
            # dataDict = {'I': I,
            #             'Q': Q,
            #             'freqs': freqs,
            #             'metadata': {'flux': flux,
            #                          'ffl_freq': qb.pars['ffl_freq'],
            #                         'ffl_LO': qb.pars['ffl_LO'],
            #                          'ffl_IF': qb.pars['ffl_IF']},
            #             }
            # save_datadict_to_fgroup(g_subcalib, 'rrSpec_ffl_drive_off', dataDict)
            inst.set_ffl_LO(qb.pars['ffl_LO'])
            
            if bOptimizeFFLMixer:
                ffl_mixer_data = qb.mixer_powers(sa, 'ffl', switch='on')
                if qb.pars['ffl_IF'] == 0.0:
                    if ffl_mixer_data['ffl_lo_leakage'] > -72:
                        qb.optimize_mixer(sa, element='ffl',cal='LO', switch='off')
                else:
                    if ffl_mixer_data['ffl_lo_leakage'] > -72:
                        qb.optimize_mixer(sa, element='ffl',cal='LO')
                    if ffl_mixer_data['ffl_im_leakage'] > -72:
                        qb.optimize_mixer(sa, element='ffl',cal='SB')
                    
               
            #mixer_data_off = qb.mixer_powers(sa, element='ffl', switch='off')
            #mixer_data_on = qb.mixer_powers(sa, element='ffl', switch='on')
            # I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=34,IF_min=200e6,IF_max=206e6,df=50e3,n_avg=4000,savedata=True)
            # dataDict = {'I': I,
            #             'Q': Q,
            #             'freqs': freqs,
            #             'metadata': {'flux': flux,
            #                          'ffl_freq': qb.pars['ffl_freq'],
            #                          'ffl_LO': qb.pars['ffl_LO'],
            #                          'ffl_IF': qb.pars['ffl_IF'],
            #                          },
            #             }
            # for k in mixer_data_on.keys():
            #     dataDict['metadata'][k + '_on'] = mixer_data_on[k]
            # for k in mixer_data_off.keys():
            #     dataDict['metadata'][k + '_off'] = mixer_data_off[k]
            # save_datadict_to_fgroup(g_subcalib, 'rrSpec_ffl_drive_on', dataDict)
            if not test: 
                sweep_powers(qb,sa,element='ffl', n_avg=n_avg, amp_ffl_scale=1, data_key=data_key,**kwargs)
                
                
                
    

def measure_base_ringdown_time(qb,sa, flux_list, n_avg =4000, rr_pulse_len=1.4e3):
    qb.update_value('rr_pulse_len_in_clk', clk(rr_pulse_len))
    inst.turn_off_ffl_drive()
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{qb.device_name}\\{sDate}\\ringdown'
    
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    filename = f'ringdownBase_sweepFLux_DA={qb.pars["rr_atten"]}dB_navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)
    
    
      
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_data = hf.create_group(f'ringdown_{timestamp}')
        mixer_data = qb.mixer_powers(sa, 'rr')
        for k in mixer_data.keys():
            g_data.attrs[k] = mixer_data[k]
        g_calib = hf.create_group(f'calib_{timestamp}')
        
        base_flux = 51e-6
        dataDict = acquire_rr_spec_background(qb, base_flux=base_flux, n_avg=n_avg, IF_min=141e6, IF_max=148e6)
        I0 = dataDict['I']
        Q0 = dataDict['Q']
        freqs = dataDict['freqs']
        z = np.polyfit(freqs, np.abs(I0 + 1j*Q0), 3)
        p = np.poly1d(z)
        save_datadict_to_fgroup(g_calib, f'flux = {base_flux*1e6:.0f} uA', dataDict)
        
        for flux in tqdm(flux_list):
            inst.set_flux_bias(flux)
            I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=141e6,IF_max=148e6,df=50e3,n_avg=n_avg,savedata=True,fit=False)
            fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q) - p(freqs))
            pf.spec_plot(freqs,I,Q,attenuation=qb.pars['rr_atten'],df=0.05e6,element='resonator',fwhm=fwhm,fc=fc, flux=flux)
            dataDict = {
                'I': I,
                'Q': Q,
                'freqs': freqs,
                'metadata': {'flux': flux,
                             'rr_freq': fc},
                }
            save_datadict_to_fgroup(g_calib, f'flux = {flux*1e6:.0f} uA', dataDict)
            qb.update_value('rr_freq', fc)
            qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'])
            
            dataDict, fig = measure_ringdown_drive_off(qb,tmax=2e3, dt=16, n_avg=n_avg, amp_r_scale = 1)
            save_datadict_to_fgroup(g_data, f'flux = {flux*1e6:.0f}', dataDict)
    
    
    
def sweep_flux(qb, sa, flux_list, n_avg=4000, rr_pulse_len=1.2e3, bCalibrateRo=True, bOptimizeFFLMixer=True, test=False):
    # frequency sweep
    stepsize = 10e6
    fmin = 2.35e9
    fmax = 3.65e9
    # save data
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{qb.device_name}\\{sDate}\\ringdown'
    
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    filename = f'ringdown_sweepFLux_sweepFFLfreq_ffl_IF={qb.pars["ffl_IF"]/1e6}MHz_amp_ffl={qb.pars["amp_ffl"]}_DA={qb.pars["rr_atten"]}dB_fDA={qb.pars["ffl_atten"]}dB_navg={n_avg}'
    index = get_index_for_filename(saveDir, filename)
        
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_data = hf.create_group(f'ringdown_{timestamp}')
        mixer_data = qb.mixer_powers(sa, 'rr')
        for k in mixer_data.keys():
            g_data.attrs[k] = mixer_data[k]
        g_calib = hf.create_group(f'calib_{timestamp}')
        
        # if not bCalibrateRo:
        #     rrTableFile = 'G:\\Shared drives\\CavityCooling\\data\\diss09_6024\\20230509\\spectroscopy\\diss_spec\\rr_freq_calibrate.csv'
        #     df = pd.read_csv(rrTableFile)
        #     fluxes = np.array(df.loc[:,'flux'])
        #     rrFreqs = np.array(df.loc[:,'rr_freq'])
        # else:
        #     base_flux = 51e-6
        #     dataDict = acquire_rr_spec_background(qb, base_flux=base_flux, n_avg=n_avg, IF_min=45e6, IF_max=54e6)
        #     I0 = dataDict['I']
        #     Q0 = dataDict['Q']
        #     freqs = dataDict['freqs']
        #     z = np.polyfit(freqs, np.abs(I0 + 1j*Q0), 3)
        #     p = np.poly1d(z)
        #     save_datadict_to_fgroup(g_calib, f'flux = {base_flux*1e6:.0f} uA', dataDict)
            
        for flux in tqdm(flux_list):
            inst.set_ffl_bias(flux, step = 10e-6, lower_bound=-10e-3, upper_bound=10e-3)
            #qb.update_value('rr_pulse_len_in_clk', clk(rr_pulse_len))
           
            if bCalibrateRo:
                inst.turn_off_ffl_drive()
                I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=29,IF_min=196e6,IF_max=204e6,df=50e3,n_avg=1000,savedata=True)
                fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
                #qb.update_value('rr_freq', fc)
                pf.spec_plot(freqs,I,Q,attenuation=qb.pars['rr_atten'],df=50e3,element='resonator',fwhm=fwhm,fc=fc, flux=flux)
                dataDict = {
                    'I': I,
                    'Q': Q,
                    'freqs': freqs,
                    'metadata': {'flux': flux,
                                 'rr_freq': fc},
                    }
                save_datadict_to_fgroup(g_calib, f'flux = {flux*1e6:.0f} uA', dataDict)
                #qb.update_value('rr_freq', fc)
                #qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'])
            # else:
            #     # read rr_freq from the look up table
            #     rr_freq = rrFreqs[np.argmin(abs(fluxes-flux*1e6))] * 1e9 
            #     qb.update_value('rr_freq', rr_freq)
            #     qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'])
            #     I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=45e6,IF_max=55e6,df=0.1e6,n_avg=n_avg,savedata=True,fit=False, fc =rr_freq)
            sweep_freq_power(qb, sa, 
                             n_avg = n_avg, 
                             stepsize = stepsize,
                             fmin = fmin,
                             fmax = fmax,
                             bOptimizeFFLMixer = bOptimizeFFLMixer, 
                             saveDir=saveDir,filename=filename, index=index, timestamp=timestamp, test=test)
        
    
    
def main():
    qb = dissipator('diss08_11a',device_name='diss08_11a')
    #qb.update_value('diss_freq', 10.7e9)
    #qb.update_value('rr_freq', 5.594e9)
    #qb.update_value('rr_LO', 5.39e9)
    #qb.update_value('rr_IF', qb.pars['rr_freq'] -qb.pars['rr_LO'] )
    #qb.update_value('rr_atten', 32)
    #qb.update_value('amp_ffl', 0.3)
#    qb.add_key('diss_freq', 9.70e9)
    #qb.update_value('ffl_freq', 2.9e9)
    #qb.update_value('ffl_LO', 2.85e9)
    #qb.update_value('ffl_LO', qb.pars['ffl_freq'])
    #qb.update_value('ffl_IF', 0.0)
    qb.update_value('ffl_atten', 9)
    #qb.update_value('_atten', 9)
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    inst.set_attenuator(qb.pars['rr_atten'])
    inst.set_rr_LO(qb.pars['rr_LO'])# turn on
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    qb.update_value('ffl_IF',300e6 )
    inst.turn_qb_LO_off()
    inst.turn_off_fflqc_drive()
    n_avg = 8000
    qb.update_value('rr_pulse_len_in_clk',1700)
    bOptimizeFFLMixer = True
    bOptimizeRRMixer = False
    bCalibrateRo = False
    test = False
        
    # flux_list = np.array([-80, -58, -49.0, -42.0, -36.0, -31.0, -26.0, -22.0, -18.0, -15.0]) * 1e-6
    flux_list = [1.15e-3]
    rr_pulse_len = 4.8e3
    
    sa = inst.init_sa()
    if bOptimizeRRMixer:
        qb.play_pulses(element='rr')
        inst.set_attenuator(0)
        rr_mixer_data = qb.mixer_powers(sa, 'rr')
        if rr_mixer_data['rr_lo_leakage'] > -75:
            qb.optimize_mixer(sa, element='rr',cal='LO')
        if rr_mixer_data['rr_im_leakage'] > -75:
            qb.optimize_mixer(sa, element='rr',cal='SB')
        
    start = timeit.default_timer() 
    
    sweep_flux(qb,sa,flux_list=flux_list,n_avg=n_avg, rr_pulse_len=rr_pulse_len, bOptimizeFFLMixer=bOptimizeFFLMixer, test=test)
    
    # measure_base_ringdown_time(qb,sa, flux_list, n_avg =n_avg, rr_pulse_len=rr_pulse_len)
    
    sa_close_device(sa)
    stop = timeit.default_timer() 
    print('Time: ', stop - start)
    
    
if __name__ == "__main__":
    main()
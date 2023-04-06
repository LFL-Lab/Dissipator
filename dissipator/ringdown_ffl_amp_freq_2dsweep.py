# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 13:58:36 2023

@author: lfl
"""

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
qb.update_value('ffl_LO', 6.550e9)
qb.update_value('ffl_IF', 50e6)
qb.update_value('amp_ffl', 0.45)
inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
bOptimizeMixer = False

from instrument_init import init_sa, init_sa_by_serial_number

if bOptimizeMixer:
    # mixer optimization
    ref_H = 20
    ref_L = -30
    qb.play_pulses()
    sa = init_sa_by_serial_number(20234492)
    qb_lo_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO'],reference=ref_H,config=True,plot=True)
    qb_im_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO']-qb.pars['ffl_IF'],reference=ref_H,config=True,plot=True)
    qb_on_power = qb.get_power(sa, freq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
    
    '''Optimize FFL mixer'''
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
    #qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element = 'ffl')
    
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

def measure_ringdown_drive_on(amp_r_scale=1, 
                              amp_ffl_scale=1,
                              tmin = 0,
                              tmax = 2e3,
                              dt = 16,
                              n_avg = 1000,freq=0.
                              ):
    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)

    # inst.set_attenuator(attenuation=qb.pars['rr_atten'])
    inst.turn_on_ffl_drive()
    resettime_clk= clk(qb.pars['rr_resettime'])
    with program() as prog:
        update_frequency('rr', (qb.pars['rr_freq']-qb.pars['rr_LO'])) 
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO']+freq)) 
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
    return dataDict, fig, fitted_pars

def measure_ringdown_drive_off(amp_r_scale=1,
                               tmin =  0,
                               tmax= 2e3,
                               dt= 16,
                               n_avg= 1000,
                               ):
    inst.turn_off_ffl_drive()
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
    
    
device = 'diss08_07A'
today = datetime.today()
sDate =  today.strftime("%Y%m%d")
expName = 'Ringdown_Sweep'
saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\ringdown'
n_avg = 3000
#%% doing sweeps
if not os.path.exists(saveDir):
    Path(saveDir).mkdir(parents=True, exist_ok=True)
filename = f'{expName}_fflFreq={str(qb.pars["ffl_freq"]/1e9).replace(".","d")}GHz_DA={qb.pars["rr_atten"]}dB_navg={n_avg}'
index = get_index_for_filename(saveDir, filename)
with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
    now = datetime.now()
    timestamp = now.strftime("%H:%M:%S")
    g_on = hf.create_group(f'sweep_ffl_amp_{timestamp}')
    amp_ffl_scale_list = np.linspace(0,0.5,10)
    ffl_IF_freq_list=np.linspace(-350e6,350e6,40)
    T_ring=np.zeros(10*40)
    i=0
    for amp_ffl_scale in amp_ffl_scale_list:
        for freq in ffl_IF_freq_list:
            dataDict, fig, fitted_pars= measure_ringdown_drive_on(amp_ffl_scale=amp_ffl_scale,tmax = 1e3,dt = 8,n_avg = 3000,freq=freq)
            T_ring[i]=fitted_pars[1]
            i=i+1
            
            
plt.imshow(T_ring.reshape((10,40)),extent=[-350,350,0.5,0],aspect="auto")
plt.colorbar();plt.xlabel('frequency (MHz)')
plt.ylabel('Amplitude scaling')
plt.show()


for i in np.arange(6):
    plt.subplot(2,3,i+1)
    plt.plot(ffl_IF_freq_list,T_ring.reshape((10,40))[i,:])
    plt.title(f'amp_scale={round(amp_ffl_scale_list[i],3)}')
    
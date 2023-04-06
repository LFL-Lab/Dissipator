# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 11:40:45 2023

@author: lfl
"""
from qubit import *
import instrument_init as inst
import h5py

qb = qubit('logical')
qb.pars['ffl_LO']=6550000000
inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
inst.set_qb_LO(qb.pars['qubit_LO'])
bOptimizeMixer = False
bCalibratePi = False

#%% mixer optimization
if bOptimizeMixer:
    inst.main()
    ref_H = 20
    ref_L = -30
    qb.play_pulses()
    # qubit mixer
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
    
    # readout mixer
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
    
    # FFL mixer
    sa_close_device(sa)
    sa = init_sa_by_serial_number(20234492)
    rr_lo_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
    rr_im_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO']-qb.pars['ffl_IF'],span = 1e6,reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power
    rr_on_power = qb.get_power(sa, freq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],reference=ref_H,config=True,plot=True) # reference should be set ABOVE expected image power

    # do a coarse sweep to minimize LO leakage
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='coarse',reference = ref_H, element='ffl')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, mode='intermediate',reference = ref_H, element='ffl')
    qb.opt_mixer( sa, cal='LO',  freq_span = 1e6, reference = ref_L, mode='fine',element='ffl')
    qb.opt_mixer( sa, cal='SB', freq_span = 1e6,  mode='coarse', reference = ref_H, element='ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='intermediate',reference = ref_H, element='ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, mode='fine', reference = ref_L, element='ffl')

# rabi calibration for pi_half length
if bCalibratePi:
    print(f'pi_amp = {qb.pars["pi_amp"]}')
    qb.update_value('pi_len', 64)
    amps, I, Q, job, period = qb.power_rabi(a_min = 0.01, a_max = 1.0, da = 0.005,  n_avg = 1000,
                                            fit = True, plot = True, detuning = 0e6, check_mixers = False)
    qb.update_value('pi_amp', round(qb.pars["pi_amp"] * period/2,3))
    
    t_arr, I, Q, job, pi2width = qb.pulse_exp(exp='rabi',n_avg = 1000, tmin = 16,tmax = 640, dt = 4, amp_q_scaling = 1, fit = True, plot = True,
                                         detuning = 0e6, check_mixers=False)
    qb.update_value('pi_half_len', int(pi2width[1]/4/4)*4)
# do background ramsey

t_arr, I, Q, job, fitted_pars = qb.pulse_exp(exp='ramsey', check_mixers = False, n_avg=5000,dt=4,tmin=qb.pars['pi_half_len'],tmax=10*200+qb.pars['pi_half_len'],detuning=0)
def measure_t2_w_ffl(amp_r_scale=1,
                    amp_ffl_scale=0.3,
                    tmin = 0,
                    tmax = 1000,
                    dt = 16,
                    n_avg = 1000,
                    detuning=0,):

    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    inst.turn_on_ffl_drive()
    resettime_clk= clk(qb.pars['qubit_resettime'])
    with program() as prog:
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']) + detuning)
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO'])) 
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    play("pi_half", "qubit")
                    align("ffl", "qubit")
                    play("pi_half", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("pi_half", "qubit")
                    align("ffl", "qubit")
                    play('const'*amp(amp_ffl_scale), "ffl", duration=t)
                    align("ffl", "qubit")
                    play("pi_half", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
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
    
    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='ramsey',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='ramsey',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                 qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                             'amp_ffl_scale': amp_ffl_scale,
                             'tmin': 16,
                             'tmax': 2e3,
                             'dt': 16,
                             'n_avg': n_avg,},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }
    return dataDict, fig


def measure_t2_ffl_off(amp_r_scale=1,
                    amp_ffl_scale=0.3,
                    tmin = 0,
                    tmax = 1000,
                    dt = 16,
                    n_avg = 1000,
                    detuning = 0,):

    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    inst.turn_off_ffl_drive()
    resettime_clk= clk(qb.pars['qubit_resettime'])
    with program() as prog:
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']) + detuning)
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    play("pi_half", "qubit")
                    align("ffl", "qubit")
                    play("pi_half", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("pi_half", "qubit")
                    wait(t, "qubit")
                    play("pi_half", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
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
    
    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='ramsey',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='ramsey',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                 qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                             'amp_ffl_scale': amp_ffl_scale,
                             'tmin': 16,
                             'tmax': 2e3,
                             'dt': 16,
                             'n_avg': n_avg,},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }
    return dataDict, fig

device = 'diss08_07A'
today = datetime.today()
sDate =  today.strftime("%Y%m%d")
saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\T2wFFL'
n_avg = 10000
if not os.path.exists(saveDir):
    Path(saveDir).mkdir(parents=True, exist_ok=True)
filename = f'T2wFFL_fflFreq={str(qb.pars["ffl_freq"]/1e9).replace(".","d")}GHz_DA={qb.pars["rr_atten"]}dB_navg={n_avg}'
index = get_index_for_filename(saveDir, filename)
with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
    now = datetime.now()
    timestamp = now.strftime("%H:%M:%S")
    g_on = hf.create_group(f'sweep_ffl_amp_{timestamp}')
    amp_ffl_scale_list = np.linspace(0.126, 0.25, 10)
    for amp_ffl_scale in amp_ffl_scale_list:
        dataDict, fig = measure_t2_w_ffl(amp_ffl_scale=amp_ffl_scale,dt = 4,tmax=1e3, n_avg=n_avg, detuning=0)
        save_datadict_to_fgroup(g_on, f'ffl amp = {amp_ffl_scale:.3f}', dataDict)

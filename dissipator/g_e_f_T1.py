# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:37:31 2023

@author: lfl
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 18:50:10 2023

@author: lfl
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 11:40:45 2023

@author: lfl
"""
from qubit import *
import instrument_init as inst
import h5py


qb = qubit('diss07a')

bOptimizeMixer = False
bCalibratePi = False

#%% mixer optimization
if bOptimizeMixer:
    inst.init_sa()
    ref_H = -10
    ref_L = -60
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
    amps, I, Q, job, period = qb.power_rabi(a_min = 0.01, a_max = 2, da = 0.005,  n_avg = 2000,fit = True, plot = True, detuning = 0e6, check_mixers = False)
    qb.update_value('pi_amp', round(qb.pars["pi_amp"] * period/2,3))
    
    t_arr, I, Q, job, pi2width = qb.pulse_exp(exp='rabi',n_avg = 2000, tmin = 16,tmax = 640, dt = 4, amp_q_scaling = 1, fit = True, plot = True,detuning = 0e6, check_mixers=False)
    qb.update_value('pi_half_len', int(pi2width[1]/4/4)*4)



def optimize_ffl_mixer(qb, sa, element):
    #qb_lo_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO'],reference=ref_H,config=True,plot=True)
    #qb_im_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO']-qb.pars['ffl_IF'],reference=ref_H,config=True,plot=True)
    #qb_on_power = qb.get_power(sa, freq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
    ref_H=0
    ref_L=-30
    '''Optimize FFL mixer'''
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H,amp_q=1, mode='coarse', element=element)
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H,amp_q=1, mode='intermediate', element=element)
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L,amp_q=1, mode='fine', element = element)
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H,amp_q=0.6, mode='coarse', element = element)
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H,amp_q=0.5, mode='intermediate', element =element)
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L,amp_q=0.5, mode='fine',element =element)

def measure_leakage_w_ffl(amp_r_scale=1,
                    amp_ffl_scale=0.3,
                    tmin = 0,
                    tmax = 5e3,
                    dt = 64,
                    n_avg = 1000,):

    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    inst.turn_on_ffl_drive()
    resettime_clk= clk(qb.pars['qubit_resettime']*2)
    with program() as prog:
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']))
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO'])) 
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    #play('const'*amp(amp_ffl_scale), "ffl", duration=1e3)
                    #align("ffl", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play('const'*amp(amp_ffl_scale), "ffl", duration=t)
                    align("ffl", "rr")
                    #wait(t, "rr")
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
    
    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='T1',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='T1',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                 qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                             'amp_ffl_scale': amp_ffl_scale,
                             'tmin': tmin,
                             'tmax': tmax,
                             'dt': dt,
                             'n_avg': n_avg,},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }
    return dataDict, fig


def measure_t1_ffl_off(amp_r_scale=1,
                    amp_ffl_scale=0.3,
                    tmin = 0,
                    tmax = 15e3,
                    dt = 64,
                    n_avg = 1000,):

    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    inst.turn_off_ffl_drive()
    resettime_clk= clk(qb.pars['qubit_resettime'])
    with program() as prog:
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']))
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    play("pi", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("pi", "qubit")
                    wait(t, "qubit")
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
    
    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='T1',dt=t_arr[-1]*1e-6/len(t_arr))
    fig = pf.plot_data(t_arr,ydata,sequence='T1',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                 qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=1)
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                             'amp_ffl_scale': amp_ffl_scale,
                             'tmin': tmin,
                             'tmax': tmax,
                             'dt': dt,
                             'n_avg': n_avg,},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }
    return dataDict, fig

def main():
    device = 'diss08_07A'
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    expName = 'T1wFFL'
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\{expName}'
    n_avg = 10000
    qb.update_value('rr_freq', qb.pars['rr_freq'])
    qb.update_value('rr_LO', qb.pars['rr_LO'])
    qb.update_value('rr_IF', qb.pars['rr_freq'] -qb.pars['rr_LO'] )
    qb.update_value('rr_atten', 32)
    qb.update_value('qubit_freq', 4.6383e9)
    qb.update_value('qubit_LO', 4.59e9)
    qb.update_value('qubit_IF', qb.pars['qubit_freq'] - qb.pars['qubit_LO'])
    qb.update_value('diss_freq', qb.pars['diss_freq'])
    qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['qubit_freq'])
    qb.update_value('ffl_LO', 3.841e9)
    qb.update_value('ffl_IF', qb.pars['ffl_freq'] - qb.pars['ffl_LO'])
    qb.add_key('ffl_atten', 14)
    inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
    inst.set_qb_LO(qb.pars['qubit_LO'])
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    #step_size = 100e6
    #ffl_freq_list = [qb.pars["ffl_freq"] + step_size * n for n in range(-3,3)]
    ffl_freq_list= [qb.pars['ffl_freq']]
    for ffl_freq in ffl_freq_list:
        qb.update_value('ffl_freq', ffl_freq)
        qb.update_value('ffl_IF', qb.pars['ffl_freq'] - qb.pars['ffl_LO'])
        filename = f'{expName}_fflFreq={str(qb.pars["ffl_freq"]/1e9).replace(".","d")}GHz_DA={qb.pars["rr_atten"]}dB_navg={n_avg}'
        index = get_index_for_filename(saveDir, filename)
        start = timeit.default_timer()   
        with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
            now = datetime.now()
            timestamp = now.strftime("%H:%M:%S")
            g_on = hf.create_group(f'sweep_ffl_amp_{timestamp}')
            dataDict, fig = measure_t1_ffl_off(dt=64, tmax=10e3, n_avg=n_avg)
            save_datadict_to_fgroup(g_on, f'ffl amp = off', dataDict)
            amp_ffl_scale_list = np.linspace(0., 1, 20)
            for amp_ffl_scale in amp_ffl_scale_list:
                if amp_ffl_scale<0.45:
                    dataDict, fig = measure_t1_w_ffl(amp_ffl_scale=amp_ffl_scale,dt = 64,tmax=6e3, n_avg=n_avg)
                elif 0.45<amp_ffl_scale<0.7:
                    dataDict, fig = measure_t1_w_ffl(amp_ffl_scale=amp_ffl_scale,dt = 4,tmax=2e3, n_avg=n_avg)
                else:
                    dataDict, fig = measure_t1_w_ffl(amp_ffl_scale=amp_ffl_scale,dt = 4,tmax=1e3, n_avg=n_avg)  
                save_datadict_to_fgroup(g_on, f'ffl amp = {amp_ffl_scale:.3f}', dataDict)
        stop = timeit.default_timer()
        print('Time: ', stop - start)  
    
if __name__ == '__main__':
    main()
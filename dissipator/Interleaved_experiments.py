# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 13:29:41 2023

@author: lfl
"""
## Interleaved cxavity cooling and cavity reset experiments for paper
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
import numpy as np
import interleaved_plot as ipf


def interleaved_cavity_reset(qb,
                     amp_r_scale=1,
                    amp_ffl_scale=1,
                    tmin = 0,
                    tmax = 5e3,
                    dt = 64,
                    n_avg = 1000,flux=0, detuning=0, ffl_len=250):
    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    resettime_clk= clk(qb.pars['qubit_resettime'])
    rr_IF = qb.pars['rr_IF']
    ffl_IF = qb.pars['ffl_IF']
    qb_IF=qb.pars['qubit_IF']+detuning
    with program() as prog:
        update_frequency('rr', rr_IF) 
        update_frequency('ffl', ffl_IF) 
        update_frequency('qubit' ,qb_IF)
        n, t = qb.declare_vars([int, int])
        I_b, Q_b, I_c, Q_c, I_cf,Q_cf = qb.declare_vars([fixed,fixed, fixed, fixed, fixed, fixed])
        I_bst, Q_bst, I_cst, Q_cst, I_cfst,Q_cfst, n_stream = qb.declare_streams(stream_num=7)
        
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_each_(t, t_arr):
                #play bare echo
                #play("readout"*amp(0.), "rr")
                #align('ffl','rr')
                #play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
                #play('gaussian'*amp(0.), "ffl", duration=clk(ffl_len))
                #wait(clk(80))
                #align('ffl','qubit')
                #wait(resettime_clk)
                play("pi_half", "qubit")
                with if_(t>0):
                    wait(t)
                play("pi", "qubit")
                with if_(t>0):
                    wait(t, "qubit")
                play("pi_half", "qubit")
                align("qubit","rr")
                measure("readout", "rr", None, *qb.res_demod(I_b, Q_b))
                save(I_b, I_bst)
                save(Q_b, Q_bst)
                wait(resettime_clk)
                align()
                
                #play echo with cavity drive
                
                play("readout"*amp(amp_r_scale), "rr")
                align('ffl','rr')
                
                play('gaussian'*amp(0.), "ffl", duration=clk(ffl_len))
                wait(clk(80))
                align('ffl','qubit')
                play("pi_half", "qubit")
                with if_(t>0):
                    wait(t)
                play("pi", "qubit")
                with if_(t>0):
                    wait(t, "qubit")
                play("pi_half", "qubit")
                align("qubit","rr")
                measure("readout", "rr", None, *qb.res_demod(I_c, Q_c))
                save(I_c, I_cst)
                save(Q_c, Q_cst)
                wait(resettime_clk)
                align()
                
                #play echo with cavity drive and ffl
                play("readout"*amp(amp_r_scale), "rr")
                align('ffl','rr')
                play('gaussian'*amp(amp_ffl_scale), "ffl", duration=clk(ffl_len))
                wait(clk(80))
                align('ffl','qubit')
                play("pi_half", "qubit")
                with if_(t>0):
                    wait(t)
                play("pi", "qubit")
                with if_(t>0):
                    wait(t, "qubit")
                play("pi_half", "qubit")
                align("qubit","rr")
                measure("readout", "rr", None, *qb.res_demod(I_cf, Q_cf))
                save(I_cf, I_cfst)
                save(Q_cf, Q_cfst)
                wait(resettime_clk)
              
                
                
                
                
        with stream_processing():
            I_bst.buffer(len(t_arr)).average().save("I_b")
            Q_bst.buffer(len(t_arr)).average().save("Q_b")
            I_cst.buffer(len(t_arr)).average().save("I_c")
            Q_cst.buffer(len(t_arr)).average().save("Q_c")
            I_cfst.buffer(len(t_arr)).average().save("I_cf")
            Q_cfst.buffer(len(t_arr)).average().save("Q_cf")
            n_stream.save('n')
    return prog


def interleaved_cavity_cooling(qb,
                     amp_r_scale=1,
                    amp_ffl_scale=1,
                    tmin = 0,
                    tmax = 5e3,
                    dt = 64,
                    n_avg = 1000,flux=0, detuning=0, ffl_len=250):
    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    resettime_clk= clk(qb.pars['qubit_resettime'])
    rr_IF = qb.pars['rr_IF']
    ffl_IF = qb.pars['ffl_IF']
    qb_IF=qb.pars['qubit_IF']+detuning
    with program() as prog:
        update_frequency('rr', rr_IF) 
        update_frequency('ffl', ffl_IF) 
        update_frequency('qubit' ,qb_IF)
        n, t = qb.declare_vars([int, int])
        I_b, Q_b, I_c, Q_c, I_cf,Q_cf = qb.declare_vars([fixed,fixed, fixed, fixed, fixed, fixed])
        I_bst, Q_bst, I_cst, Q_cst, I_cfst,Q_cfst, n_stream = qb.declare_streams(stream_num=7)
        
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_each_(t, t_arr):
                #play bare echo
                play('const'*amp(0.), "ffl", duration=2*t+100)
                wait(40,"qubit")
                play("pi_half", "qubit")
                align('rr','qubit')
                play("readout"*amp(0.), "rr", duration=2*t+20)
                with if_(t>0):
                    wait(t,"qubit")
                play("pi", "qubit")
                with if_(t>0):
                    wait(t,"qubit")
                play("pi_half", "qubit")
                align("ffl","rr")
                measure("readout", "rr", None, *qb.res_demod(I_b, Q_b))
                save(I_b, I_bst)
                save(Q_b, Q_bst)
                wait(resettime_clk)
                align()
                
                #play echo with cavity drive
                play('const'*amp(0.), "ffl", duration=2*t+100)
                wait(40,"qubit")
                play("pi_half", "qubit")
                align('qubit','rr')
                play("readout"*amp(amp_r_scale), "rr", duration=2*t+20)
                with if_(t>0):
                    wait(t,"qubit")
                play("pi", "qubit")
                with if_(t>0):
                    wait(t,"qubit")
                play("pi_half", "qubit")
                #align("qubit","rr")
                align("ffl","rr")
                measure("readout", "rr", None, *qb.res_demod(I_c, Q_c))
                save(I_c, I_cst)
                save(Q_c, Q_cst)
                wait(resettime_clk)
                align()
                
                
                #play echo with cavity drive and ffl
                play('const'*amp(amp_ffl_scale), "ffl", duration=2*t+100)
                wait(40,"qubit")
                play("pi_half", "qubit")
                align('rr','qubit')
                play("readout"*amp(amp_r_scale), "rr", duration=2*t+20)
                with if_(t>0):
                    wait(t,"qubit")
                play("pi", "qubit")
                with if_(t>0):
                    wait(t,"qubit")
                play("pi_half", "qubit")
                #align("qubit","rr")
                align("ffl","rr")
                measure("readout", "rr", None, *qb.res_demod(I_cf, Q_cf))
                save(I_cf, I_cfst)
                save(Q_cf, Q_cfst)
                wait(resettime_clk)
                
                
        with stream_processing():
            I_bst.buffer(len(t_arr)).average().save("I_b")
            Q_bst.buffer(len(t_arr)).average().save("Q_b")
            I_cst.buffer(len(t_arr)).average().save("I_c")
            Q_cst.buffer(len(t_arr)).average().save("Q_c")
            I_cfst.buffer(len(t_arr)).average().save("I_cf")
            Q_cfst.buffer(len(t_arr)).average().save("Q_cf")
            n_stream.save('n')
    return prog


def interleaved_exp(qb ,sa = 0,
                  exp='interleaved_cavity_reset',
                  check_mixers=False,
                  n_avg = 2000,
                  tmin = 16,         # minimum pulse duration in nanoseconds
                  tmax = 10e3,    # maximum pulse duration in nanoseconds
                  dt = 500,        # step of sweep in nanoseconds
                  amp_q_scaling = 1,
                  fit = True,
                  plot = 'single',
                  detuning = 0e6, play_init_pi=True, simulate=True, ffl_len=1e3, amp_ffl_scale=0, scrambling_amp=1, with_ffl=True, amp_r_scale=1, flux=0):
    """

    Args:
        exp (str):   What experiment to run. Options are: 'rabi','ramsey','T1','echo'.
        n_avg (TYPE, optional): DESCRIPTION. Defaults to 2000.
        t0 (TYPE, optional): minimum pulse duration in nanoseconds. Defaults to 0.
        tf (TYPE, optional): maximum pulse duration in nanoseconds. Defaults to 10e3.
        dt (TYPE, optional): sequence step size in nanoseconds. Defaults to 500.
        amp_q_scaling (TYPE, optional): DESCRIPTION. Defaults to 1.
        plot (TYPE, optional): Whether to plot and fit the data. Defaults to True.
        detuning (float): detuning from fL0-fIF in Hz.
        resettime (TYPE, optional): waiting time between experiments. Defaults to 400e3.

    Returns:
        times (TYPE): DESCRIPTION.
        I (TYPE): DESCRIPTION.
        Q (TYPE): DESCRIPTION.

    """

    try:
        list_of_files = glob.glob(f'{saveDir}\{exp}\*.csv')
        latest_file = max(list_of_files, key=os.path.getctime)
        iteration = int(latest_file[-7:-4].lstrip('0')) + 1
    except:
        iteration = 1

    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)

    inst.set_attenuator(attenuation=qb.pars['rr_atten'])
    inst.set_ffl_attenuator(attenuation=qb.pars['ffl_atten'])

    # qb.update_value('qubit_LO', qb.pars['qubit_LO'])
    # qb.update_value('rr_LO', qb.pars['rr_LO'])

    # qb_IF = qb.pars['qubit_freq']-qb.pars['qubit_LO']
    # qb.update_value('qubit_IF', (qb.pars['qubit_freq']-qb.pars['qubit_LO']) + detuning)
    #qb.check_mix_cal(sa, amp_q = amp_q_scaling,check = check_mixers, threshold = -55)

    if exp=='interleaved_cavity_reset':
        prog = interleaved_cavity_reset(qb, amp_r_scale=amp_r_scale, amp_ffl_scale=amp_ffl_scale,tmin = 4*tmin,tmax = 4*tmax,dt = 4*dt,n_avg = n_avg,flux=flux, detuning=detuning, ffl_len=ffl_len ) 
    elif exp=='interleaved_cavity_cooling':
        prog =  interleaved_cavity_cooling(qb, amp_r_scale=amp_r_scale, amp_ffl_scale=amp_ffl_scale,tmin = 4*tmin,tmax = 4*tmax,dt = 4*dt,n_avg = n_avg,flux=flux, detuning=detuning, ffl_len=ffl_len )
    
    
    if simulate:
        qmm = QuantumMachinesManager(host=host, port=port)
        job = qmm.simulate(config=qb.config, program=prog, simulate=SimulationConfig(duration=10000))
        job.get_simulated_samples().con1.plot()
        
    datadict, job = qb.get_results(jobtype = prog, result_names = ["I_b", "Q_b", "I_c", "Q_c", "I_cf", "Q_cf"], n_total=n_avg, notify = False)
    
    

    # qb_power = qb.get_power(sa,freq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],reference=0,amp_q = amp_q_scaling, span=1e6,config=True,output=False)

    # t_arr = np.array(datadict["t"])/1e3 # times in microseconds
    
    t_arr = np.array(t_arr)*4/1e3 * 2   ##Echo therefore multiplied by 2

    
    I_b = np.array(datadict["I_b"])
    Q_b = np.array(datadict["Q_b"])
    I_c = np.array(datadict["I_c"])
    Q_c = np.array(datadict["Q_c"])
    I_cf = np.array(datadict["I_cf"])
    Q_cf = np.array(datadict["Q_cf"])
    
    # ydata_b = np.abs(I_b+1j*Q_b)
    # ydata_c = np.abs(I_c+1j*Q_c)
    # ydata_cf = np.abs(I_cf+1j*Q_cf)
    ydata_b = I_b
    ydata_c = I_c
    ydata_cf = I_cf
    if fit:
        fitted_pars_b, error_b = pf.fit_data(t_arr,ydata_b,sequence=exp,dt=t_arr[-1]*1e-6/len(t_arr))
        fitted_pars_c, error_c = pf.fit_data(t_arr,ydata_c,sequence=exp,dt=t_arr[-1]*1e-6/len(t_arr))
        fitted_pars_cf, error_cf = pf.fit_data(t_arr,ydata_cf,sequence=exp,dt=t_arr[-1]*1e-6/len(t_arr))
        if plot == 'sep':
            fig = pf.plot_data(t_arr,ydata_b,sequence=exp,fitted_pars=fitted_pars_b,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                     qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF']+detuning,fflDriveFreq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],iteration=iteration,amp_ffl_scale=amp_ffl_scale, amp=amp_r_scale,flux=flux, error=error_b, ffl_len=ffl_len, rr_atten=qb.pars['rr_atten'], ffl_atten=qb.pars['ffl_atten'])
            fig = pf.plot_data(t_arr,ydata_c,sequence=exp,fitted_pars=fitted_pars_c,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                     qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF']+detuning,fflDriveFreq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],iteration=iteration,amp_ffl_scale=amp_ffl_scale, amp=amp_r_scale,flux=flux, error=error_c, ffl_len=ffl_len, rr_atten=qb.pars['rr_atten'], ffl_atten=qb.pars['ffl_atten'])
            fig = pf.plot_data(t_arr,ydata_cf,sequence=exp,fitted_pars=fitted_pars_cf,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
                     qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF']+detuning,fflDriveFreq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],iteration=iteration,amp_ffl_scale=amp_ffl_scale, amp=amp_r_scale,flux=flux, error=error_cf, ffl_len=ffl_len, rr_atten=qb.pars['rr_atten'], ffl_atten=qb.pars['ffl_atten'])
        elif plot=='single':
            fig = ipf.plot_data_inter(t_arr,ydata_b,ydata_c,ydata_cf,seq=exp, qubitDriveFreq=qb.pars['qubit_freq']+detuning,fflDriveFreq=qb.pars['ffl_freq'],
                                          nAverages=n_avg,
                                          stepSize=dt,fitted_pars_c=fitted_pars_c,fitted_pars_cf=fitted_pars_cf,fitted_pars_b=fitted_pars_b,
                                          savefig=True, amp=amp_r_scale, ffl_atten=qb.pars['ffl_atten'], rr_atten=qb.pars['rr_atten'], flux=flux, amp_ffl_scale=amp_ffl_scale, error_b=error_b,error_c=error_c,error_cf=error_cf, ffl_len=ffl_len)
            


    #print(error)


    # if exp == 'rabi':
        # qb.update_value('pi_half_len',4*(clk(fitted_pars[1]/4)))

    # qb.update_value('qubit_IF',qb_IF)

    exp_dict = {'date/time':     datetime.now(),
               'nAverages': n_avg,
                     'Tmax': tmax,
                     'dt':   dt,
                     'pi2': qb.pars['pi_half_len'],
                     'A_d':     amp_q_scaling,
                     'w_d':  qb.pars['qubit_LO']+qb.pars['qubit_IF']+detuning,
                     'w_LO': qb.pars['qubit_LO'],
            'wait_period':  qb.pars['qubit_resettime'],
            'ffl_freq': qb.pars['ffl_freq'],
            'pulse_len': qb.pars['rr_pulse_len_in_clk'],
            'amp_ffl_scale':  amp_ffl_scale,
            'amp_r_scale':  amp_r_scale,
            'flux': flux
            }

    # save data
    dataPath = f'{saveDir}\\{exp}'
    if not os.path.exists(dataPath):
        Path(dataPath).mkdir(parents=True, exist_ok=True)
    with open(f"{dataPath}\data_{iteration:03d}.csv","w") as datafile:
        writer = csv.writer(datafile)
        writer.writerow(exp_dict.keys())
        writer.writerow(exp_dict.values())
        writer.writerow(t_arr)
        writer.writerow(I_b)
        writer.writerow(Q_b)
        writer.writerow(I_c)
        writer.writerow(Q_c)
        writer.writerow(I_cf)
        writer.writerow(Q_cf)
    fig.savefig(f"{dataPath}\data_{iteration:03d}.png")
    
    dataDict = {
        'I_b': I_b,
        'Q_b': Q_b,
        'I_c': I_c,
        'Q_c': Q_c,
        'I_cf': I_cf,
        'Q_cf': Q_cf,
        'fitted_pars_b': fitted_pars_b,
        'fitted_pars_c': fitted_pars_c,
        'fitted_pars_cf': fitted_pars_cf,
        'errorb': error_b,
        'errorc': error_c,
        'errorcf': error_cf,
        't_arr': t_arr,
        'metadata': {'flux': flux},
                    
        }
    
    return dataDict

def main():
    
    qb = dissipator('diss08_11a',device_name='diss08_11a')

    qb.update_value('ffl_atten', 9)
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    inst.set_attenuator(qb.pars['rr_atten'])
    inst.set_rr_LO(qb.pars['rr_LO'])# turn on
    inst.set_qb_LO(qb.pars['qubit_LO'])
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    n_avg = 12000
    amp_ffl_scale=0.
    amp_r_scale=0.9
    flux=inst.get_ffl_bias()*1e3
    
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{qb.device_name}\\{sDate}\\interleaved_exp'
    
    
    expName='interleaved_cavity_cooling'
    
    
    filename = f'{expName}_sweepffllen_flux={flux*1e3}mA_fflFreq={(qb.pars["ffl_freq"])/1e9:.2f}GHz_DA={qb.pars["rr_atten"]}dB_fDA={qb.pars["ffl_atten"]}dB_amp_ffl_scale={amp_ffl_scale}_navg={n_avg}amp_r_scale={amp_r_scale}'
    
    
    
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    
    index = get_index_for_filename(saveDir, filename)
    
    with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
        now = datetime.now()
        timestamp = now.strftime("%H:%M:%S")
        g_on = hf.create_group(f'sweep_ffl_atten_cav_amp', track_order=True)
    
        #ffl_len_list = np.flip(np.concatenate((4*np.round(np.linspace(4, 120, 15)), 4*np.round(np.linspace(120,750, 15)))))
        #ffl_len_list = np.array([96])
        #ffl_len_list=np.flip(4*np.round(np.linspace(24,660,22)))
        #atten_list=[10,14,18,22,26]
        atten_list=[9,9,9,9,9,9,9,9,9,9,9,9,9]
        #for cav_amp in [0.]:
        #g_subdata = g_on.create_group(f'amp_r_scale={amp_r_scale:.3f}', track_order=True)
        if expName=='interleaved_cavity_cooling':
            for j, atten in enumerate(atten_list):
                qb.update_value('ffl_atten', atten)
                inst.set_ffl_attenuator(qb.pars['ffl_atten'])
                dataDict= qb.pulse_exp(exp = 'cavity-cooling', n_avg = n_avg, tmin =40, tmax = 12e3, dt = 64, fit=True, check_mixers=False, simulate=False, with_ffl=True, amp_ffl_scale=1., amp_r_scale=0., detuning=0, flux=inst.get_ffl_bias()*1e3) 
                save_datadict_to_fgroup(g_on, f'ffl iter = {j:.1f}', dataDict)
        else:
            for j, ffl_len in enumerate(ffl_len_list):
                
                #dataDict= interleaved_exp(qb, exp=expName, n_avg = n_avg,tmin = 40,tmax = 12e3, dt = 64,amp_q_scaling = 1,simulate=True, ffl_len=1300, amp_ffl_scale=1., scrambling_amp=1, with_ffl=True, amp_r_scale=cav_amp, flux=flux)
                #save_datadict_to_fgroup(g_subdata, f'ffl len = {ffl_len:.1f}', dataDict)
                if ffl_len<1000:
                    dataDict=qb.pulse_exp(exp = 'cavity-reset', n_avg = n_avg, tmin =20, tmax = 3e3, dt = 4, fit=True, check_mixers=False, simulate=False, ffl_len=ffl_len, with_ffl=True, amp_ffl_scale=0., amp_r_scale=0.9, detuning=0, flux=inst.get_ffl_bias()*1e3)
                else:
                    dataDict=qb.pulse_exp(exp = 'cavity-reset', n_avg = n_avg, tmin =40, tmax = 12e3, dt = 64, fit=True, check_mixers=False, simulate=False, ffl_len=ffl_len, with_ffl=True, amp_ffl_scale=0., amp_r_scale=0.9, detuning=0, flux=inst.get_ffl_bias()*1e3)
                save_datadict_to_fgroup(g_on, f'ffl len = {ffl_len:.1f}', dataDict)
                
        
    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    main()

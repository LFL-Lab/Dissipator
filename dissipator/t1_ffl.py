# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 18:50:10 2023

@author: lfl
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 11:40:45 2023

@author: lfl

qubit reset as a function of frequency and amplitude
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




#%% mixer optimization
def optimize_mixer(sa, qb, element='rr', cal='LO'):
    ref_H = -10
    ref_L = -50
    qb.play_pulses(element)
    if element == 'rr':
        inst.set_attenuator(0)
        inst.get_attenuation()
    qb_lo_leakage = qb.get_power(sa, freq=qb.pars[f'{element}_LO'],reference=ref_H,config=True,plot=True)
    qb_im_leakage = qb.get_power(sa, freq=qb.pars[f'{element}_LO']-qb.pars[f'{element}_IF'],reference=ref_H,config=True,plot=True)
    qb_on_power = qb.get_power(sa, freq=qb.pars[f'{element}_LO']+qb.pars[f'{element}_IF'],reference=ref_H, config=True,plot=True)# reference should be set ABOVE expected image power
    
    '''Optimize FFL mixer'''
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='coarse', element = element)
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='coarse', element = element)
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='intermediate', element = element)
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='intermediate', element = element)
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_L, mode='fine', element = element)
    qb.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_L, mode='fine', element = element)
    
    
    qb.write_pars()


# rabi calibration for pi_half length
def calibrate_pi_pulse(qb):
    print(f'pi_amp = {qb.pars["pi_amp"]}')
    qb.update_value('pi_len', 64)
    amps, I, Q, job, period = qb.power_rabi(a_min = 0.01, a_max = 2, da = 0.005,  n_avg = 2000,fit = True, plot = True, detuning = 0e6, check_mixers = False)
    qb.update_value('pi_amp', round(qb.pars["pi_amp"] * period/2,3))
    
    t_arr, I, Q, job, pi2width = qb.pulse_exp(exp='rabi',n_avg = 10000, tmin = 16,tmax = 640, dt = 4, amp_q_scaling = 1, fit = True, plot = True,detuning = 0e6, check_mixers=False)
    qb.update_value('pi_half_len', int(pi2width[1]/4/4)*4)




def measure_t1_w_ffl(qb,
                     amp_r_scale=1,
                    amp_ffl_scale=1,
                    amp_fflqc_scale=0.3,
                    tmin = 0,
                    tmax = 5e3,
                    dt = 64,
                    n_avg = 1000,flux=0):

    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    #inst.turn_on_ffl_drive()
    
    resettime_clk= clk(qb.pars['qubit_resettime'])
    with program() as prog:
        update_frequency('rr', qb.pars['rr_IF'])
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']))
        update_frequency('ffl', qb.pars['ffl_IF']) 
        #update_frequency('fflqc', qb.pars['fflqc_IF']) 
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_each_(t,t_arr):
                with if_(t==0):
                    play("pi", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
                    wait(resettime_clk)
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("pi", "qubit")
                    #align("qubit", "rr")
                    wait(clk(60),'ffl')
                    play('gaussian'*amp(amp_ffl_scale), "ffl", duration=t)
                    wait(clk(30))
                    align("ffl","rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
                    wait(resettime_clk)
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
                  qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1, amp_ffl_scale=amp_ffl_scale, ffl_atten=qb.pars['ffl_atten'],rr_atten=qb.pars['rr_atten'], flux=inst.get_ffl_bias()*1e3, error=error)
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
    return dataDict, fig, fitted_pars
    #return prog

def measure_t1_w_ffl_qd_df(qb,
                     amp_r_scale=1,
                    amp_ffl_scale=1.0,
                    amp_fflqc_scale=0.3,
                    tmin = 0,
                    tmax = 5e3,
                    dt = 64,
                    n_avg = 500,
                    ffl_len = 2e3):

    tmin = clk(tmin)
    tmax = clk(tmax)
    dt = clk(dt)
    t_arr = np.arange(tmin, tmax + dt/2, dt, dtype = int)
    inst.turn_on_fflqc_drive()
    
    resettime_clk= clk(qb.pars['qubit_resettime'])
    with program() as prog:
        update_frequency('rr', qb.pars['rr_IF'])
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']))
        update_frequency('ffl', qb.pars['ffl_IF']) 
        #update_frequency('fflqc', qb.pars['fflqc_IF']) 
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            # with for_each_(t,t_arr):
                
                # with if_(t==0):
                #     play("pi", "qubit")
                #     align("qubit", "rr")
                #     measure("readout", "rr", None, *qb.res_demod(I, Q))
                #     wait(resettime_clk)
                #     save(I, I_stream)
                #     save(Q, Q_stream)
                    
                # with else_():
            play("pi", "qubit")
            align("ffl", "qubit")
            #align("fflqc", "qubit")
            #align('qubit','rr')
            play('gaussian'*amp(amp_ffl_scale), "ffl", duration=ffl_len)
            #play('gaussian'*amp(amp_fflqc_scale), "fflqc", duration=ffl_len)
            #align("qubit","fflqc")
            #align("fflqc","rr")
            align("ffl", "rr")
            # wait(t, "rr")
            measure("readout", "rr", None, *qb.res_demod(I, Q))
            wait(resettime_clk)
            save(I, I_stream)
            save(Q, Q_stream)

        with stream_processing():
            # I_stream.buffer(len(t_arr)).average().save("I")
            # Q_stream.buffer(len(t_arr)).average().save("Q")
            I_stream.average().save("I")
            Q_stream.average().save("Q")
            n_stream.save('n')
    
    
    
    datadict, job = qb.get_results(jobtype = prog, result_names = ["I", "Q"], n_total=n_avg, notify = False)
    t_arr = np.array(t_arr)*4/1e3
    I = np.array(datadict["I"])
    Q = np.array(datadict["Q"])
    ydata = np.abs(I+1j*Q)
    
    # fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='T1',dt=t_arr[-1]*1e-6/len(t_arr))
    # fig = pf.plot_data(t_arr,ydata,sequence='T1',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
    #               qubitDriveFreq=qb.pars['fflqc_freq'],qb_power = -8,iteration=1, amp=amp_fflqc_scale, ffl_atten=qb.pars['ffl_atten'],rr_atten=qb.pars['rr_atten'])
    
    #print(ydata)
    dataDict = {'metadata': {'amp_r_scale': amp_r_scale,
                              'amp_ffl_scale': amp_fflqc_scale,
                              'tmin': tmin,
                              'tmax': tmax,
                              'dt': dt,
                              'n_avg': n_avg,},
                'time': t_arr,
                'I': I,
                'Q': Q,
                
        }
    return ydata
    # return prog

# host='10.71.0.56'
# port='9510'
# qmm = QuantumMachinesManager(host=host, port=port)
# job = qmm.simulate(config=qb.config, program = measure_t1_w_ffl(qb, n_avg=2, tmin=100), simulate=SimulationConfig(duration=500))
# job.get_simulated_samples().con1.plot()



def measure_t1_ffl_off(qb,
                       amp_r_scale=1,
                    amp_ffl_scale=0.0,
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
        update_frequency('rr', qb.pars['rr_IF'])
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']))
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_each_(t,t_arr):
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
                  qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1, amp_ffl_scale=amp_ffl_scale, ffl_atten=qb.pars['ffl_atten'],rr_atten=qb.pars['rr_atten'], flux=inst.get_ffl_bias()*1e3, error=error)
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

def ffl_punchout_amp_sweep(qb):
    #qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['qubit_freq'])
    qb.update_value('ffl_IF', 150e6)
    qb.update_value('ffl_LO',qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    stepsize = 2
    atten_list = np.arange(30,0,-2)
    amp_list= np.linspace(0,1,15)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    data = []
    inst.turn_on_ffl_drive()  
    for amp in amp_list:
        inst.set_ffl_attenuator(16)
        I, Q, freqs, job = resonator_spec_wffl(qb,f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=20e6,IF_max=70e6,amp_ffl_scale=amp,df=0.1e6,n_avg=2000,savedata=True)
        data.append(np.abs(I + 1j * Q))
        
        
        
    #im=ax.imshow(data, aspect='auto',origin='lower',extent=(freqs[0]/1e9, freqs[-1]/1e9,atten_list[0], atten_list[-1]),
            #interpolation=None, cmap='RdBu')
    im=ax.imshow(data, aspect='auto',origin='lower',extent=(freqs[0]/1e9, freqs[-1]/1e9, amp_list[0], amp_list[-1]),
            interpolation=None, cmap='RdBu')
    #ax.legend()
    fig.colorbar(im)
    plt.show()
    return fig

def ffl_punchout(qb, sweep_mode='atten', ffl_on=True, stepsize = 2) :
    #qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['qubit_freq'])
    # qb.update_value('ffl_IF', 350e6)
    # qb.update_value('ffl_LO',qb.pars['ffl_freq'] - qb.pars['ffl_IF'])
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    data = []
    inst.turn_on_ffl_drive()  
    if sweep_mode == 'atten':
        atten_list = np.arange(30,0,-stepsize)
        for atten in atten_list:
            inst.set_ffl_attenuator(atten)
            if ffl_on:
                I, Q, freqs, job = resonator_spec_wffl(qb,f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=20e6,IF_max=70e6,amp_ffl_scale=0.,df=0.1e6,n_avg=2000,savedata=True)
            else:
                I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=20e6,IF_max=70e6,df=0.1e6,n_avg=2000,savedata=True)
            data.append(np.abs(I + 1j * Q))
        im=ax.imshow(data, aspect='auto',origin='lower',extent=(freqs[0]/1e9, freqs[-1]/1e9,atten_list[0], atten_list[-1]),
                    interpolation=None, cmap='RdBu')
    elif sweep_mode =='amp'        :
        amp_list= np.linspace(0,1,10)
        for amp in amp_list:
            inst.set_ffl_attenuator(9)
            I, Q, freqs, job = resonator_spec_wffl(qb,f_LO=qb.pars['rr_LO'],atten=qb.pars['rr_atten'],IF_min=35e6,IF_max=55e6,amp_ffl_scale=amp,df=0.1e6,n_avg=5000,savedata=True)
            data.append(np.abs(I + 1j * Q))
        im=ax.imshow(data, aspect='auto',origin='lower',extent=(freqs[0]/1e9, freqs[-1]/1e9,amp_list[0], amp_list[-1]),
                interpolation=None, cmap='RdBu')
    
    ax.set_title(f'FFL f_d={round(qb.pars["ffl_freq"]/1e9,5)}GHz')
    fig.colorbar(im)
    plt.show()
    return fig

def resonator_spec_wffl(qb, IF_min = 0.1e6,
                       f_LO = 7e9,
                       IF_max = 400e6,
                       df = 0.1e6,
                       atten = 10,
                       n_avg = 500,
                       amp_ffl_scale=1,
                       res_ringdown_time = int(5e3),
                       port_type = 'notch',
                       fit=True,
                       savedata=True):

    try:
        list_of_files = glob.glob(r'D:\weak_measurements\spectroscopy\resonator_spec\*.csv')
        latest_file = max(list_of_files, key=os.path.getctime)
        iteration = int(latest_file[-7:-4].lstrip('0')) + 1
    except:
        iteration = 1

    freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
    # freqs_list = freqs.tolist()
    # set attenuation and change rr_LO freq
    inst.set_attenuator(attenuation=atten)
    qb.update_value('rr_LO', value = f_LO)
    inst.set_rr_LO(qb.pars['rr_LO'])

    ### QUA code ###
    with program() as rr_spec:

        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)
        n_stream = declare_stream()
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO'])) 

        with for_(n, 0, n < n_avg, n + 1):

            # with for_each_(f, freqs_list):
            with for_(f, IF_min, f < IF_max + df/2, f + df):
                
                update_frequency("rr", f)
                wait(res_ringdown_time, "rr")
                align('rr','ffl')
                play('const'*amp(amp_ffl_scale), "ffl", duration=clk(3e3))
                #align("ffl", "rr")
                measure("readout", "rr", None,*qb.res_demod(I, Q))
                save(I, I_st)
                save(Q, Q_st)

            save(n,n_stream)

        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')
            n_stream.save('n')

    datadict,job = qb.get_results(rr_spec,result_names=["I","Q","n"],showprogress=False)

    I = np.array(datadict["I"])
    Q = np.array(datadict["Q"])
    freq_arr = np.array(freqs + qb.pars['rr_LO'])

    if fit:
        fc,fwhm = pf.fit_res(freq_arr,np.abs(I+1j*Q))
        pf.spec_plot(freq_arr,I,Q,attenuation=atten,df=df,iteration=iteration,element='resonator',fwhm=fwhm,fc=fc)
        print(f'Resonant Frequency: {fc*1e-9:.5f} GHz\nFWHM = {fwhm*1e-6} MHz\nkappa = {2*np.pi*fwhm*1e-6:.3f} MHz')

    exp_dict = {'date/time':    datetime.now(),
               'nAverages': n_avg,
                     'w_LO': qb.pars['rr_LO'],
            'wait_period':  res_ringdown_time,
            }
    if savedata:
        # save data
        dataPath = f'{saveDir}\spectroscopy\\resonator_spec'
        if not os.path.exists(dataPath):
            Path(dataPath).mkdir(parents=True, exist_ok=True)
        with open(f"{dataPath}\data_{iteration:03d}.csv","w") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(exp_dict.keys())
            writer.writerow(exp_dict.values())
            writer.writerow(freqs)
            writer.writerow(I)
            writer.writerow(Q)

    return I, Q, freqs+qb.pars['rr_LO'], job;




def main():
    
    qb = dissipator('diss08_11a',device_name='diss08_11a')
    #qb.update_value('diss_freq', 10.7e9)
    #qb.update_value('rr_freq', 5.5928e9)
    #qb.update_value('rr_LO', 5.39e9)
    #qb.update_value('rr_IF', qb.pars['rr_freq'] -qb.pars['rr_LO'] )
    #qb.update_value('rr_atten', 32)
    #qb.update_value('amp_ffl', 0.4)
    #    qb.add_key('diss_freq', 9.70e9)
    #qb.update_value('ffl_freq', 2.216e9)
    #qb.update_value('ffl_LO', 2.85e9)
    #qb.update_value('ffl_LO', 5.06e9)
    qb.update_value('ffl_IF', 0)
    #qb.update_value('fflqc_LO', 2.217e9)
    #qb.update_value('fflqc_IF', 0)
    qb.update_value('ffl_atten', 4)
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    inst.set_attenuator(qb.pars['rr_atten'])
    #qb.update_value('qubit_freq', 3.377e9)
    #qb.update_value('qubit_LO', 3.1e9)
    #qb.update_value('qubit_IF',qb.pars['qubit_freq'] - qb.pars['qubit_LO'] )
    inst.set_rr_LO(qb.pars['rr_LO'])# turn on
    #inst.set_ffl_LO_off(qb.pars['ffl_LO'])
    inst.set_qb_LO(qb.pars['qubit_LO'])
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    n_avg = 2000
    bOptimizeFFLMixer = True
    bOptimizeRRMixer = False
    bCalibrateRo = True
    test = False
    sa = inst.init_sa()     
    # flux_list = np.array([-80, -58, -49.0, -42.0, -36.0, -31.0, -26.0, -22.0, -18.0, -15.0]) * 1e-6
    flux_list = np.array([15])*1e-6
    # fmin=5.05e9
    # fmax=5.05e9
    # stepsize=50e6
    # num_of_points = int((fmax - fmin)/stepsize) + 1
    #freqs_list = [fmin + n * stepsize for n in range(num_of_points)]
    freqs_list= np.linspace(4.9,5.3,50)*1e9
    start = timeit.default_timer() 
    t1vals= np.zeros((len(freqs_list),5))
    for flux in flux_list:
        rr_pulse_len = 4e3
        #sa_close_device(sa)
        
        if bOptimizeRRMixer:
            qb.play_pulses(element='rr')
            inst.set_attenuator(0)
            rr_mixer_data = qb.mixer_powers(sa, 'rr')
            if rr_mixer_data['rr_lo_leakage'] > -75:
                qb.optimize_mixer(sa, element='rr',cal='LO')
            if rr_mixer_data['rr_im_leakage'] > -75:
                qb.optimize_mixer(sa, element='rr',cal='SB')
            
        #inst.set_ffl_bias(flux, step = 5e-6, lower_bound=-1e-3, upper_bound=5e-3)
        
        if bCalibrateRo:
            inst.turn_off_ffl_drive()
            inst.turn_off_fflqc_drive()
            I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=30,IF_min=200e6,IF_max=206e6,df=50e3,n_avg=1000,savedata=True)
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
            qb.update_value('rr_freq', fc)
            qb.update_value('rr_IF', qb.pars['rr_freq'] - qb.pars['rr_LO'])
            #qb.update_value('fflqc_LO', qb.pars['rr_freq']-qb.pars['qubit_freq'])
            #qb.update_value('fflqc_IF', 0)
            #qb.update_value('fflqc_freq', qb.pars['rr_freq']-qb.pars['qubit_freq'])
            #inst.set_fflqc_LO(qb.pars['fflqc_LO'])
        for i, ffl_freq in enumerate(freqs_list):
            qb.update_value('ffl_LO', ffl_freq)
            qb.update_value('ffl_freq', ffl_freq)
            inst.set_ffl_LO(ffl_freq)
            inst.turn_on_ffl_drive()
            qb.update_value('ffl_IF', 0)
            #qb.update_value('ffl_freq', ffl_lo+qb.pars['ffl_IF'])
            if bOptimizeFFLMixer:
                ffl_mixer_data = qb.mixer_powers(sa, 'ffl', switch='off')
                if ffl_mixer_data['ffl_lo_leakage'] > -67:
                    qb.optimize_mixer(sa, element='ffl',cal='LO', switch='off')
                    qb.mixer_powers(sa, 'ffl')
                #optimize_mixer(sa, qb, element='ffl',cal='SB')
            # for ffl_IF in IF_list:
            #     qb.update_value('ffl_IF', ffl_IF)
            #     qb.update_value('ffl_freq', ffl_lo + ffl_IF)
       
            #     if bOptimizeFFLMixer:
                    
            #         optimize_mixer(sa, qb, element='ffl',cal='SB')
            #     # sweep amplitude
            
            rrLen = qb.pars['rr_pulse_len_in_clk']
            expName= 't1ffl'
            filename = f'{expName}_sweepPowers_flux={flux*1e3}mA_fflFreq={(qb.pars["ffl_freq"])/1e9:.2f}GHz_DA={qb.pars["rr_atten"]}dB_fDA={qb.pars["ffl_atten"]}dB_rrLen={rrLen}clks_navg={n_avg}'
            index = get_index_for_filename(saveDir, filename)
            inst.set_attenuator(qb.pars['rr_atten'])
            with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
                now = datetime.now()
                timestamp = now.strftime("%H:%M:%S")
                g_on = hf.create_group(f'sweep_ffl_amp_{timestamp}')
                    
                dataDict, fig = measure_t1_ffl_off(qb, dt=1024, tmin=96,tmax=60e3, n_avg=n_avg)
                save_datadict_to_fgroup(g_on, f'ffl amp = off', dataDict)
            
                amp_ffl_scale_list = np.linspace(0., 1, 3)
                for j, amp_ffl_scale in enumerate(amp_ffl_scale_list):
                    dataDict, fig, fits = measure_t1_w_ffl(qb,amp_ffl_scale=amp_ffl_scale,tmin=96, dt = 1024,tmax=60e3, n_avg=n_avg)
                    print(fits)
                    save_datadict_to_fgroup(g_on, f'ffl amp = {amp_ffl_scale:.3f}', dataDict)
                    t1vals[i,j]=fits[1]
    sa_close_device(sa)
    #qb.update_value('rr_pulse_len_in_clk', 2500)
    stop = timeit.default_timer()
    print('Time: ', stop - start)  
    
if __name__ == '__main__':
    main()

# fig, ax = plt.subplots(figsize=(6,6))
# plot=ax.imshow(t1vals, interpolation='nearest', extent=[0,1,4.7,5.4], vmin=2, vmax=24, aspect="auto", origin="lower")
# #ax.set_yticklabels([250, 500,1000,2000, 3000, 5000])
# plt.xlabel("ffl amp")
# plt.ylabel("freq")
# plt.colorbar(plot)
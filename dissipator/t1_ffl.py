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
from qubit import *
import instrument_init as inst
import h5py
import timeit



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
                    play("pi", "qubit")
                    align("qubit", "rr")
                    measure("readout", "rr", None, *qb.res_demod(I, Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("pi", "qubit")
                    align("ffl", "qubit")
                    play('const'*amp(amp_ffl_scale), "ffl", duration=t)
                    align("ffl", "rr")
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
                 qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1)
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


def measure_t1_ffl_off(qb,
                       amp_r_scale=1,
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
                 qubitDriveFreq=qb.pars['ffl_freq'],qb_power = -8,iteration=1)
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
                play('const'*amp(amp_ffl_scale), "ffl", duration=clk(4e3))
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
    
    qb = qubit('diss07a')

    bOptimizeRRMixer = False
    bOptimizeFFLMixer = True
    bCalibratePi = False
    
    device = 'diss08_07A'
    today = datetime.today()
    sDate =  today.strftime("%Y%m%d")
    expName = 'T1wFFL'
    saveDir = f'G:\\Shared drives\\CavityCooling\data\\{device}\\{sDate}\\{expName}'
    n_avg = 10000
    
    

    qb.update_value('ffl_freq', qb.pars['diss_freq'] - qb.pars['qubit_freq'])
    qb.update_value('ffl_IF', 150e6)
    qb.update_value('ffl_LO', qb.pars['ffl_freq'] -qb.pars['ffl_IF'])
    # qb.update_value('ffl_IF', qb.pars['ffl_freq'] - qb.pars['ffl_LO'])
    qb.add_key('ffl_atten', 20)
  
    inst.set_rr_LO(qb.pars['rr_LO'])
    inst.set_qb_LO(qb.pars['qubit_LO'])
    inst.set_ffl_attenuator(qb.pars['ffl_atten'])
    if not os.path.exists(saveDir):
        Path(saveDir).mkdir(parents=True, exist_ok=True)
    
    start = timeit.default_timer() 
    # mixer calibration
    sa = inst.init_sa()
    
    if bOptimizeRRMixer:
        qb.play_pulses(element='rr')
        inst.set_attenuator(0)
        optimize_mixer(sa, qb, element='rr',cal='LO')
        optimize_mixer(sa, qb, element='rr',cal='SB')
        
    # IF_min = 50e6
    # IF_max = 250e6
    # stepsize = 50e6
    # IF_list = [int(qb.pars['ffl_IF'] + stepsize * n) for n in range(-4,3)]
        
    #step_size = 100e6
    #ffl_freq_list = [qb.pars["ffl_freq"] + step_size * n for n in range(-3,3)]
    ffl_freq_list= qb.pars['ffl_freq']
    ffl_freq_list = [3.3417e9, 3.4417e9,3.5417e9,3.2417e9,3.1417e9,3.0417e9,2.9417e9,2.8417e9,2.7417e9,2.6417e9]
    for ffl_freq in ffl_freq_list:
        qb.update_value('ffl_freq', ffl_freq)
        qb.update_value('ffl_LO', qb.pars['ffl_freq']-qb.pars['ffl_IF'])
        inst.set_ffl_LO(qb.pars['ffl_LO'])
        #qb.update_value('ffl_freq', ffl_lo+qb.pars['ffl_IF'])
        if bOptimizeFFLMixer:
            optimize_mixer(sa, qb, element='ffl',cal='LO')
            optimize_mixer(sa, qb, element='ffl',cal='SB')
        # for ffl_IF in IF_list:
        #     qb.update_value('ffl_IF', ffl_IF)
        #     qb.update_value('ffl_freq', ffl_lo + ffl_IF)
   
        #     if bOptimizeFFLMixer:
                
        #         optimize_mixer(sa, qb, element='ffl',cal='SB')
        #     # sweep amplitude
        
        rrLen = qb.pars['rr_pulse_len_in_clk']
        filename = f'{expName}_sweepPowers_flux=10uA_fflFreq={(qb.pars["ffl_freq"])/1e9:.2f}GHz_DA={qb.pars["rr_atten"]}dB_fDA={qb.pars["ffl_atten"]}dB_rrLen={rrLen}clks_navg={n_avg}'
        index = get_index_for_filename(saveDir, filename)
        with h5py.File(f'{saveDir}\\{filename}_{index}.h5','w') as hf:
            now = datetime.now()
            timestamp = now.strftime("%H:%M:%S")
            g_on = hf.create_group(f'sweep_ffl_amp_{timestamp}')
                
            dataDict, fig = measure_t1_ffl_off(qb, dt=64, tmax=10e3, n_avg=n_avg)
            save_datadict_to_fgroup(g_on, f'ffl amp = off', dataDict)
        
            amp_ffl_scale_list = np.linspace(0., 0.99, 7)
            for amp_ffl_scale in amp_ffl_scale_list:
                dataDict, fig = measure_t1_w_ffl(qb,amp_ffl_scale=amp_ffl_scale,dt = 8,tmax=4e3, n_avg=n_avg)
                save_datadict_to_fgroup(g_on, f'ffl amp = {amp_ffl_scale:.3f}', dataDict)
    sa_close_device(sa)
    stop = timeit.default_timer()
    print('Time: ', stop - start)  
    
if __name__ == '__main__':
    main()
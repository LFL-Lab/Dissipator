import timeit
from qubit import *
import instrument_init as inst
import plot_functions as pf
import h5py
from datetime import datetime
import os
from qubit import qubit

qb = qubit('logical')
qb.update_value('ffl_LO', 6.550e9)
qb.update_value('ffl_IF', 50e6)
qb.update_value('amp_ffl', 0.45)
inst.set_ffl_LO(qb.pars['ffl_LO']) # turn on
bOptimizeMixer = False
if bOptimizeMixer:
    # mixer optimization
    inst.main()
    ref_H = 20
    ref_L = -30
    qb.play_pulses()
    qb_lo_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO'],reference=ref_H,config=True,plot=True)
    qb_im_leakage = qb.get_power(sa, freq=qb.pars['ffl_LO']-qb.pars['ffl_IF'],reference=ref_H,config=True,plot=True)
    qb_on_power = qb.get_power(sa, freq=qb.pars['ffl_LO']+qb.pars['ffl_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
    
    '''Optimize mixer'''
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
    qb.opt_mixer(sa, cal='LO', freq_span = 1e6, reference = ref_L, mode='fine', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='coarse', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_H, mode='intermediate', element = 'ffl')
    qb.opt_mixer(sa, cal='SB', freq_span = 1e6, reference = ref_L, mode='fine',element = 'ffl')


def measure_t1_drive_on(amp_r_scale=1, 
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

    # inst.set_attenuator(attenuation=qb.pars['rr_atten'])
    inst.turn_on_ffl_drive()
    resettime_clk= clk(qb.pars['rr_resettime'])
    with program() as prog:
        update_frequency('rr', (qb.pars['rr_freq']-qb.pars['rr_LO'])) 
        
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    play("readout"*amp(amp_r_scale), "rr")
                    play('const'*amp(amp_ffl_scale), "ffl", duration=4*qb.pars["rr_pulse_len_in_clk"])
                    align("ffl","rr")
                    measure("readout_diss", "rr", None,*qb.res_demod_diss(I,Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("readout"*amp(amp_r_scale), "rr")
                    play('const'*amp(amp_ffl_scale), "ffl", duration=4*qb.pars["rr_pulse_len_in_clk"])
                    align("rr", "ffl")
                    wait(t, "rr")
                    measure("readout_diss", "rr", None,*qb.res_demod_diss(I,Q))
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

    fitted_pars, error = pf.fit_data(t_arr,ydata,sequence='dissT1',dt=t_arr[-1]*1e-6/len(t_arr))
    pf.plot_data(t_arr,ydata,sequence='T1_drive=ON',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
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
    return dataDict

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
        
        n, t, I, Q = qb.declare_vars([int, int, fixed, fixed])

        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)

        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with for_(*from_array(t,t_arr)):
                with if_(t==0):
                    play("readout"*amp(amp_r_scale), "rr")
                    measure("readout_diss", "rr", None,*qb.res_demod_diss(I,Q))
                    wait(resettime_clk, "rr")
                    save(I, I_stream)
                    save(Q, Q_stream)
                with else_():
                    play("readout"*amp(amp_r_scale), "rr")
                    wait(t, "rr")
                    measure("readout_diss", "rr", None,*qb.res_demod_diss(I,Q))
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
    pf.plot_data(t_arr,ydata,sequence='ringdown_drive=off',fitted_pars=fitted_pars,nAverages=n_avg, pi2Width=qb.pars['pi_half_len'],
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
    
    return dataDict

def save_datadict_to_fgroup(f, name, datadict):
    subgroup = f.create_group(name)
    dset_i = subgroup.create_dataset('I', data=datadict['I'])
    dset_q = subgroup.create_dataset('Q', data=datadict['Q'])
    dset_t = subgroup.create_dataset('t', data=datadict['time'])
    for key in datadict['metadata'].keys():
        subgroup.attrs[key] = datadict['metadata'][key]
    print(f'write dataset to {name}')
    
 
def main():
    rr_atten = 20
    # do res spec
    # I, Q, freqs, job = qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=rr_atten,IF_min=30e6,IF_max=60e6,df=0.1e6,n_avg=1000,savedata=True)
    # fc,fwhm = pf.fit_res(freqs,np.abs(I+1j*Q))
    qb.update_value('rr_freq', 5.6807e9)
    # dataDict = measure_t1_drive_off(amp_r_scale=1, tmax=4e3, dt=64, n_avg=5000)
    dataDict = measure_t1_drive_on(amp_r_scale=1, amp_ffl_scale=0.3, n_avg=1000000, dt = 4, tmax=400)
    

if __name__ == "__main__":
    main()
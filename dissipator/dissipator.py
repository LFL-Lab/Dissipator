# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 00:12:50 2023

@author: lfl
"""
from qubit import *
import h5py
from datetime import datetime
import os
from instrument_init import init_sa, init_sa_by_serial_number

class dissipator(qubit):
    
    def __init__(self, qb,device_name):
        qubit.__init__(self, qb)
        self.device_name = device_name
        
    def optimize_mixer(self,sa, element='rr', cal='LO'):
        ref_H = 0
        ref_L = -45
        self.play_pulses(element)
        if element == 'rr':
            inst.set_attenuator(0)
            inst.get_attenuation()
        qb_lo_leakage = self.get_power(sa, freq=self.pars[f'{element}_LO'],reference=ref_H,config=True,plot=True)
        qb_im_leakage = self.get_power(sa, freq=self.pars[f'{element}_LO']-self.pars[f'{element}_IF'],reference=ref_H,config=True,plot=True)
        qb_on_power = self.get_power(sa, freq=self.pars[f'{element}_LO']+self.pars[f'{element}_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
        
        if qb_lo_leakage > -75: 
            self.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='coarse', element = element)
        if qb_im_leakage > -75:
            self.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='intermediate', element = element)
            self.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_H, mode='intermediate', element = element)
        # self.opt_mixer(sa, cal=cal, freq_span = 1e6, reference = ref_L, mode='fine', element = element)
        
        self.write_pars()
        
    def mixer_powers(self, sa, element='rr'):
        ref_H = 0
        ref_L = -45
        qb_lo_leakage = self.get_power(sa, freq=self.pars[f'{element}_LO'],reference=ref_H,config=True,plot=True)
        qb_im_leakage = self.get_power(sa, freq=self.pars[f'{element}_LO']-self.pars[f'{element}_IF'],reference=ref_H,config=True,plot=True)
        qb_on_power = self.get_power(sa, freq=self.pars[f'{element}_LO']+self.pars[f'{element}_IF'],reference=ref_H, config=True,plot=True) # reference should be set ABOVE expected image power
        
        return {f'{element}_lo_leakage': qb_lo_leakage, 
                f'{element}_im_leakage': qb_im_leakage,
                f'{element}_on_power': qb_on_power,}
    
    def ffl_spec(self,
                 element = 'ffl',
                 f_LO = 7e9,
                 IF_min = 0.1e6,
                 IF_max = 400e6,
                 df = 0.1e6,
                 amp_ffl_scale = 1,
                 amp_r_scale = 1,
                 n_avg = 2000,
                 savedata = True,
                 check_mixers = False,
                 plot = True):
        # TODO: modify for direct spec
        # save data
        today = datetime.today()
        sDate =  today.strftime("%Y%m%d")
        saveDir = f'G:\\Shared drives\\CavityCooling\data\\{self.device_name}\\{sDate}'
        
        dataPath = f'{saveDir}\\spectroscopy\\{element}_spec'
        filename = 'data'
        iteration = get_index_for_filename(dataPath, filename, file_format='csv')
            
        freqs = np.arange(IF_min, IF_max + df/2, df, dtype=int)
        # freqs_list = freqs.tolist()
        # set attenuation and change rr_LO freq
        inst.set_attenuator(attenuation=self.pars['rr_atten'])
        inst.set_ffl_attenuator(self.pars['ffl_atten'])
        inst.set_rr_LO(self.pars['rr_LO'])
        self.update_value('ffl_LO', f_LO)
        inst.set_ffl_LO(self.pars['ffl_LO'])
        self.update_value('ffl_IF', int(IF_max/2))
        if check_mixers:
            sa = inst.init_sa()
            self.play_pulses(element='ffl')
            self.optimize_mixer(sa, element='ffl',cal='LO')
            self.optimize_mixer(sa, element='ffl',cal='SB')
            sa_close_device(sa)
        resettime_clk = clk(self.pars['rr_resettime'])
        ### QUA code ###
        with program() as ffl_spec:
    
            n = declare(int)
            I = declare(fixed)
            I_st = declare_stream()
            Q = declare(fixed)
            Q_st = declare_stream()
            f = declare(int)
            n_stream = declare_stream()
            
            I_b, Q_b, I_tot,Q_tot = self.declare_vars([fixed,fixed, fixed, fixed])
            
            update_frequency("rr", self.pars['rr_freq'] - self.pars['rr_LO'])
            with for_(n, 0, n < n_avg, n + 1):
                with for_(f, IF_min, f < IF_max + df/2, f + df):
                    update_frequency("ffl", f)
                    
                    measure("readout", "rr", None, *self.res_demod(I_b, Q_b))
                    wait(resettime_clk, "rr")
                    align("rr", "ffl")
                    play('const'*amp(amp_ffl_scale), "ffl", duration=self.pars['rr_pulse_len_in_clk'])
                    measure("readout", "rr", None,*self.res_demod(I,Q))
                    wait(resettime_clk, "rr")
                    
                    assign(I_tot, I-I_b)
                    assign(Q_tot, Q-Q_b)
                    save(I_tot, I_st)
                    save(Q_tot, Q_st)
    
                save(n,n_stream)
    
            with stream_processing():
                I_st.buffer(len(freqs)).average().save('I')
                Q_st.buffer(len(freqs)).average().save('Q')
                n_stream.save('n')
    
        datadict,job = self.get_results(ffl_spec,result_names=["I","Q","n"],n_total=n_avg,progress_key = "n", showprogress=True)
    
        I = np.array(datadict["I"])
        Q = np.array(datadict["Q"])
        freq_arr = np.array(freqs + self.pars['ffl_LO'])
    
        if plot:
            pf.spec_plot(freq_arr,I,Q,attenuation=self.pars['ffl_atten'],df=df,iteration=iteration,element='ffl')
    
        exp_dict = {'date/time':    datetime.now(),
                   'nAverages': n_avg,
                         'ffl_LO': self.pars['ffl_LO'],
                'wait_period': self.pars['rr_resettime'],
                'report': str(job.execution_report()),
                }
        if savedata:
            # save data
            if not os.path.exists(dataPath):
                Path(dataPath).mkdir(parents=True, exist_ok=True)
            with open(f"{dataPath}\data_{iteration:03d}.csv","w") as datafile:
                writer = csv.writer(datafile)
                writer.writerow(exp_dict.keys())
                writer.writerow(exp_dict.values())
                writer.writerow(freqs)
                writer.writerow(I)
                writer.writerow(Q)
    
        return I, Q, freq_arr, job;
    
    
   
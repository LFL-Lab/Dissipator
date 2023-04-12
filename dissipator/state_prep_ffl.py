# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 23:42:02 2023

@author: lfl

state prep w and w/o ffl
"""

bScramble = False
bFFLdrive = False
qb.update_value('qubit_reset_pulse_len', 64)
qb.

def make_sequence(bScramble=False, bResete=False):
    resetpulse_clk = clk(qb.pars['qubit_reset_pulse_len'])
    with program() as prog:
        update_frequency('qubit', (qb.pars['qubit_freq']-qb.pars['qubit_LO']))
        update_frequency('ffl', (qb.pars['ffl_freq']-qb.pars['ffl_LO'])) 
        a, n, I, Q = qb.declare_vars([fixed, int, fixed, fixed])
    
        I_stream, Q_stream, n_stream = qb.declare_streams(stream_num=3)
    
        with for_(n, 0, n < n_avg, n + 1):
            save(n, n_stream)
            with if_(bScramble):
                play("scamble", "qubit")
            with if_(bReset):
                play("reset", "ffl")
            with else_():
                wait(resetpulse_clk, "ffl")
            with for_(*from_array(a,var_arr)):
                play("pi" * amp(a), "qubit")
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



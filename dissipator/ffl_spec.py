# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:13:40 2023

@author: lfl
"""

import timeit
from dissipator import *
import instrument_init as inst
import plot_functions as pf
import h5py
from datetime import datetime
import os
from instrument_init import init_sa, init_sa_by_serial_number



def main():
    qb = dissipator('diss09', device_name='diss09_5578')
    start = timeit.default_timer()
    inst.set_ffl_LO(qb.pars['ffl_LO'])
    bOptimizerrMixer = True
    if bOptimizerrMixer:
        sa = inst.init_sa()
        qb.play_pulses(element='rr')
        qb.optimize_mixer(sa, element='rr',cal='LO')
        qb.optimize_mixer(sa, element='rr',cal='SB')
        sa_close_device(sa)
    # I, Q, freq, job = qb.ffl_spec(f_LO = qb.pars['ffl_LO'])
    I, Q, freq, job = qb.run_scan(element='ffl',lo_min = 2e9, lo_max=4e9, chunksize=400e6, amp_q_scaling=0.9)
    stop = timeit.default_timer()
    print('Time: ', stop - start)  
if __name__ == "__main__":
    main()
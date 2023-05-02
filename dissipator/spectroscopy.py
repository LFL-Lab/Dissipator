# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 16:27:20 2023

@author: lfl
"""

"""
logical qubit and dissiptor spectroscopy
"""
from dissipator import *
import instrument_init as inst
device = 'diss09_5578'
qb = dissipator('diss09', device_name=device)

qb.resonator_spec(f_LO=qb.pars['rr_LO'],atten=18,IF_min=63e6,IF_max=93e6,df=0.1e6,n_avg=20000,savedata=True)

# qb.punchout(n_avg=200,span = 30e6, atten_range = [5,30], atten_step = 2, res_freq = [5.58e9])
#  
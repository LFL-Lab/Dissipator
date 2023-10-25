# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:31:33 2023

@author: lfl
"""

from qubit import qubit
import instrument_init as inst
import plot_functions as pf
import numpy as np

qb = dissipator('diss08_11a',device_name='diss08_11a')

inst.set_rr_LO(qb.pars['rr_LO'])
inst.set_qb_LO(qb.pars['qubit_LO'])
inst.set_attenuator(attenuation=qb.pars['rr_atten'])

prog = qb.make_sequence(exp='IQblob', n_reps = 2000)

datadict, job = qb.get_results(jobtype = prog,result_names=['n', 'I','Q','Iexc','Qexc'], showprogress=True, liveplot = False)

 
plot, ax = pf.init_IQ_plot()
pf.plot_single_shot(datadict)#,axes=ax)


# optimize frequency

# optimize amplitude

# optimize intergration time
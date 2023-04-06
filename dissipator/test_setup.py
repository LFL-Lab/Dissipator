# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 12:15:31 2023

@author: lfl
"""

from qubit import qubit
'''Initialize qubit class'''
qb = qubit('logical')

'''Update important parameters'''
qb.update_value('qubit_LO', value = 4.5e9)
qb.update_value('qubit_IF', value = 50e6)
qb.update_value('rr_freq', value = 5.681e9)
qb.update_value('rr_IF', 45e6)
qb.update_value('rr_LO', value = qb.pars['rr_freq'] - qb.pars['rr_IF'])

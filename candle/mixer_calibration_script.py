# Tested on Keysight FieldFox N9917A

from qm.QuantumMachinesManager import QuantumMachinesManager
import time
import numpy as np
import scipy.optimize as opti
from mixer_calib_funcs import get_power,opt_mixer

##### Qubit mixer calibration #####
mixer1 = 'qubit'
# get LO leakage power
qb_lo_leakage = get_power(sa, freq=qb_LO,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get image leakage power
qb_im_leakage = get_power(sa, freq=qb_LO-qb_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get qubit drive power at (almost) top of fridge
qb_on_power = get_power(sa, freq=qb_LO+qb_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power

# do a coarse sweep to minimize LO leakage
opt_mixer(sa, cal='lo', mode='coarse',element=mixer1)
# do a finer sweep in desired range
opt_mixer(sa, cal='lo', mode='fine',element=mixer1)
# do a coarse sweep to minimize sideband
opt_mixer(sa, cal='sb', mode='coarse',element=mixer1)
# do a finer sweep in desired range
opt_mixer(sa, cal='sb', mode='fine',element=mixer1)

##### do IQ imbalance calibration #####
mixer2 = 'rr'
# get LO leakage power
rr_lo_leakage = get_power(sa, freq=rr_LO,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get image leakage power
rr_im_leakage = get_power(sa, freq=rr_LO-rr_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# get qubit drive power at (almost) top of fridge
rr_on_power = get_power(sa, freq=rr_LO+rr_IF,reference=0,config=True,plot=True) # reference should be set ABOVE expected image power
# do a coarse sweep to minimize LO leakage
opt_mixer(sa, cal='lo', mode='coarse',element=mixer2)
# do a finer sweep in desired range
opt_mixer( sa, cal='lo', mode='fine',element=mixer2)
# do a coarse sweep to minimize sideband
opt_mixer( sa, cal='sb', mode='coarse',element=mixer2)
# do a finer sweep in desired range
opt_mixer(sa, cal='sb', mode='fine',element=mixer2)

sa_close_device(sa)

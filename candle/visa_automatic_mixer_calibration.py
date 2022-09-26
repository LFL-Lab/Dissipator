# Tested on Keysight FieldFox N9917A

from qm.QuantumMachinesManager import QuantumMachinesManager
import time
import numpy as np
import scipy.optimize as opti
from instruments_MUN11 import *



mixer = 'qubit'
qmm = QuantumMachinesManager()
qm = qmm.open_qm(config)
calib = mixer_calibration(qm,mixer)


##### do DC offset calibration #####
# get LO leakage power
freqs, signals = calib.read_off_leakage_power(plot=True, reference=0) # reference should be set ABOVE expected image power

# get signal
calib.config_sa_sweep(calib.LO_freq, reference=0) # config the sa 

# do a coarse sweep
values, argmin = calib.brute_force_search_dc_offsets(plot=True)

# do a finer sweep in desired range
values, argmin = calib.brute_force_search_dc_offsets(Imin=-0.01, Imax=0.01, Qmin= -0.02, Qmax=0.02, num_of_points=30, plot=True)


##### do IQ imbalance calibration #####
# get image power
freqs, signals = calib.read_off_image_power(plot=True, reference=0) # reference should be set ABOVE expected image power

# get signal
calib.config_sa_sweep(calib.LO_freq-calib.IF_freq,reference=0) # config the sa 

# do a coarse sweep
values, argmin = calib.brute_force_search_imbalance(plot=True)

# do a finer sweep in desired range
values, argmin = calib.brute_force_search_imbalance(pmin=-0.07, pmax=0.0, gmin= -0.02, gmax=0.05, num_of_points=30, plot=True)





# fig, axes = plt.subplots(2,2)
# calib.sweep_and_plot(calib.LO_freq,  reference=0, ax=axes[0,0])
# leakage0 = calib.get_amp()


# Parameters for Nelder-Mead
initial_simplex = np.zeros([3, 2])
initial_simplex[0, :] = [0, 0]
initial_simplex[1, :] = [0, 0.1]
initial_simplex[2, :] = [0.1, 0]
xatol = 1e-4  # 1e-4 change in DC offset or gain/phase
fatol = 3  # dB change tolerance
maxiter = 50  # 50 iterations should be more then enough, but can be changed.

# Optimize LO leakage
start_time = time.time()
if mixer == 'qubit':
    # i0 = config["controllers"]["con1"]["analog_outputs"][1]["offset"]
    # q0 = config["controllers"]["con1"]["analog_outputs"][2]["offset"]
    i0 =-0.011151
    q0 = -0.00885
elif mixer == 'readout':
    i0 = config["controllers"]["con1"]["analog_outputs"][3]["offset"]
    q0 = config["controllers"]["con1"]["analog_outputs"][4]["offset"]
fun_leakage = lambda x: calib.get_leakage(x[0], x[1])
bnds = ((-0.1, 0.1), (-0.1, 0.1))
res_leakage = opti.minimize(
    fun_leakage,
    [i0, q0],
    method="Nelder-Mead",
    bounds = bnds,
    options={
        "xatol": xatol,
        "fatol": fatol,
        "initial_simplex": initial_simplex,
        "maxiter": maxiter,
        
    },
)

print(
    f"LO Leakage Results: Found a minimum of {int(res_leakage.fun)} dBm at I0 = {res_leakage.x[0]:.5f}, Q0 = {res_leakage.x[1]:.5f} in "
    f"{int(time.time() - start_time)} seconds --- {leakage0 - int(res_leakage.fun)} dBc"
)

calib.sweep_and_plot(calib.LO_freq, ax=axes[0,1])


# Optimize image
calib.config_sa_sweep(calib.LO_freq - calib.IF_freq,reference=-50) # config the sa 
calib.sweep_and_plot(calib.LO_freq - calib.IF_freq, ax=axes[1,0])
image0 = calib.get_amp()
start_time = time.time()
if mixer == 'qubit':
    g0 = 0.06
    p0 = -0.110
elif mixer == 'readout':
    g0 = -0.012
    p0 = -0.047
fun_image = lambda x: calib.get_image(x[0], x[1])
bnds = ((-0.15, 0.15), (-0.15, 0.15))
res_image = opti.minimize(
    fun_image,
    [0, 0],
    method="Nelder-Mead",
    bounds = bnds,
    options={
        "xatol": xatol,
        "fatol": fatol,
        "initial_simplex": initial_simplex,
        "maxiter": maxiter,
    },
)
print(
    f"Image Rejection Results: Found a minimum of {int(res_image.fun)} dBm at g0 = {res_image.x[0]:.5f}, p0 = {res_image.x[1]:.5f} in "
    f"{int(time.time() - start_time)} seconds --- {image0 - int(res_image.fun)} dBc"
)

calib.sweep_and_plot(calib.LO_freq - calib.IF_freq, ax=axes[1,1])
# calib.close()
fig.tight_layout()
fig.show()

### use the following code to check at the RF port of the IQ mixer ###
fig2, ax = plt.subplots()
istar, qstar = res_leakage.x[0], res_leakage.x[1]
gstar, pstar = res_image.x[0], res_image.x[1]
# istar, qstar =  -0.00955, 0.00356
# gstar, pstar = 0.00110, -0.06841
# calib = mixer_calibration(qm,mixer)
# set mixer at the optimal bias point
# calib.qm.set_dc_offset_by_qe(calib.element, "I", istar)
# calib.qm.set_dc_offset_by_qe(calib.element, "Q", qstar)
calib.qm.set_mixer_correction(calib.mixer,int(calib.IF_freq), int(calib.LO_freq), calib.IQ_imbalance_correction(gstar, pstar))
# calib.sweep_and_plot(calib.LO_freq,reference=-20,span=200e6, sweepBW=100e3, ax = ax) # config the sa 
fig2.show()
calib.close()

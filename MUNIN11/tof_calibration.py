from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from config_MUN11 import config
import matplotlib.pyplot as plt
import numpy as np

# communicate with the server, and create qmm API
qmm = QuantumMachinesManager()

# open a new quantum machine on the server, and create an qm API
# qm = qmm.open_qm(config)

# create a QUA program, that will later be excuted on the qm
# with program() as test:g
#     play("const", "qubit")
    
# simulate the QUA prog, and create a job API
# job = qmm.simulate(config, test, SimulationConfig(1000))
# samples = job.get_simulated_samples()
# samples.con1.plot()

# ###################
# # The QUA program #
# ###################
with program() as tof_cal:
        
    n = declare(int)
    
    adc_st = declare_stream(adc_trace=True)
    
    # measure("readout", "rr", adc_st)
    

    with for_(n, 0, n < 100, n + 1):
        reset_phase("rr")
        measure("readout", "rr", adc_st)
        wait(50000, "rr")
        
    with stream_processing():
        adc_st.input1().save("adc1")
        adc_st.input2().save("adc2")

# ######################################
# # Open Communication with the Server #
# ######################################
# qmm = QuantumMachinesManager()

# ####################
# # Simulate Program #
# ####################
# simulation_config = SimulationConfig(
#                     duration=2000,
#                     simulation_interface=LoopbackInterface([("con1", 3, "con1", 1), ("con1", 4, "con1", 2)]))

# job = qmm.simulate(config, tof_cal, simulation_config)
qm = qmm.open_qm(config)
job = qm.execute(tof_cal)
res_handles = job.result_handles
res_handles.wait_for_all_values()
adc1 = res_handles.get("adc1").fetch_all()
adc2 = res_handles.get("adc2").fetch_all()

# # plt.figure()
# # plt.title('time-of-flight calibration sequence')
# # job.get_simulated_samples().con1.plot()
# # plt.xlabel("time [ns]")
# # plt.ylabel("DAC [V]")

plt.figure()
plt.title('time-of-flight calibration analysis')
plt.plot(adc1)
plt.plot(adc2)
plt.legend(["adc1", "adc2"])
print(np.mean(adc1)/4096)
print(np.mean(adc2)/4096)
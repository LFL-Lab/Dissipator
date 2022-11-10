from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from config import config
import matplotlib.pyplot as plt
import numpy as np
from VISAdrivers import LO845m as LO
import json

# communicate with the server, and create qmm API
qmm = QuantumMachinesManager()

with open('pars.json', 'r') as openfile:
    pars_dict = json.load(openfile)

with program() as tof_cal:

    n = declare(int)
    adc_st = declare_stream(adc_trace=True)
    update_frequency('rr',10e6)
    with for_(n, 0, n < 1000, n + 1):

        reset_phase("rr")
        measure("readout", "rr", adc_st)
        wait(50000, "rr")

    with stream_processing():
        adc_st.input1().average().save("adc1")
        adc_st.input2().average().save("adc2")

qm = qmm.open_qm(config)
job = qm.execute(tof_cal)
res_handles = job.result_handles
res_handles.wait_for_all_values()
adc1 = res_handles.get("adc1").fetch_all()
adc2 = res_handles.get("adc2").fetch_all()

plt.figure()
plt.title('time-of-flight calibration analysis')
plt.plot(adc1)
plt.plot(adc2)
plt.legend(["adc1", "adc2"])
offset1 = np.mean(adc1)/4096
offset2 = np.mean(adc2)/4096
print(f'Input 1 Offset: {offset1*1e3} mV')
print(f'Input 2 Offset: {offset2*1e3} mV')
pars_dict['analog_input_offsets'][0] = pars_dict['analog_input_offsets'][0] - offset1
pars_dict['analog_input_offsets'][1] = pars_dict['analog_input_offsets'][1] - offset2
with open("pars.json", "w") as outfile:
    json.dump(pars_dict, outfile)

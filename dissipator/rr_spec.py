from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from qualang_tools.addons.InteractivePlotLib import InteractivePlotLib 
from config_DISS01 import *
import matplotlib.pyplot as plt
import numpy as np
from analysis_functions import *
from file_utility import *
from datetime import datetime
from utilities import *

###################
# The QUA program #
###################
# add save

def run_rr_spec(qm, f_min = -35e6, f_max = -10e6, n_avg=10, df =0.05e6, plot=True):
    now = datetime.now()
    path = f"G:\Shared drives\LFL\Projects\PhotonVac\DISS08_07A\{now.strftime('20%y%m')}\\rr_spec\\"
    create_data_directory(path)
    Qplt = InteractivePlotLib(path)
    
    freqs = np.arange(f_min, f_max + df/2, df, dtype=int)
    freqs_list = freqs.tolist()
    readLO = rr_LO
    

    with program() as rr_spec:
    
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)
    
        # I_test = declare(fixed)
        # I_test_st = declare_stream()
        # Q_test = declare(fixed)
        # Q_test_st = declare_stream()    
    
        with for_(n, 0, n < n_avg, n + 1):

            with for_each_(f,freqs_list):
                update_frequency("rr", f)
                wait(25000, "rr")

                measure("readout", "rr", None,
                        dual_demod.full("integW_cos", "out1", "integW_sin", "out2", I),
                        dual_demod.full("integW_minus_sin", "out1", "integW_cos", "out2", Q))
                save(I, I_st)
                save(Q, Q_st)

    
        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')



    job = qm.execute(rr_spec)
    res_handles = job.result_handles
    res_handles.wait_for_all_values()
    I2 = res_handles.get("I").fetch_all()
    Q2 = res_handles.get("Q").fetch_all()
        
    if plot:

        Qplt.figure(1)
        # plt.clf()
        plt.plot(freqs+readLO,(np.sqrt(I2**2 + Q2**2)))
        plt.title("dual demod")
        plt.xlabel('freq (GHz)')
        plt.ylabel('Log Mag (mV)')
        plt.title(f'Readout Pulse Amplitude = {config["waveforms"]["ro_wf1"]["sample"]}V')
        plt.show()
        
    return freqs+readLO, I2, Q2
if __name__ == "__main__":
    qmm = QuantumMachinesManager()
    qmm.close_all_quantum_machines()
    qm = qmm.open_qm(config)
    f, I,Q = run_rr_spec(qm, f_min=0.1e6,f_max=90e6,n_avg=2000)


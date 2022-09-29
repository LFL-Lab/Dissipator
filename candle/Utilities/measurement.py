from tqdm import tqdm

from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager

from slacker import sendslack
from Utilities.data import *
   


"""
getIQ(
      config,       :: configuration dictionary (usually imported from a file named config.py),
                       which is interpretted by QUA to determine settings on OPX 
      jobtype):     :: QUA jobtype (e.g. qubit_spec, rr_spec)
"""
def getIQ(config, jobtype, showprogress = False, n_avg = 1000, notify = False):
    
    
    # Open Communication with the Server
    qmm = QuantumMachinesManager()

    # execute the job
    qm = qmm.open_qm(config)
    job = qm.execute(jobtype)
    res_handles = job.result_handles

    # get handles for I and Q
    I_handle = res_handles.get("I")
    Q_handle = res_handles.get("Q")
    I_handle.wait_for_values(1)
    Q_handle.wait_for_values(1)
    
    # display progress bar and send slack notification
    if showprogress:
        
        # get handle for iteration number n, and initialize counter
        n_handle = res_handles.get("n")
        n0 = 0 
        
        # create progressbar using tqdm
        with tqdm(total = n_avg) as progress_bar:
            
            while(I_handle.is_processing()):
                
                # retrieve iteration value n
                n = n_handle.fetch_all()
                
                # send slack notification after 100 iterations
                if notify and n0 <= 100 and n > 100:
                    
                    later = gettimedelta(progress_bar)
                    sendslack(message = f":zap: Measurement running. The measurement will be complete at approximately " + later + ".")
                
                # update progressbar with increase in n, and reset counter
                Δn = n - n0
    
                if Δn > 0:
                    progress_bar.update(Δn)
                    n0 = n
                    
               
                    
    # retrieve I and Q values      
    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    
    # close quantum machine
    qmm.close_all_quantum_machines()
    
    return I, Q, job;




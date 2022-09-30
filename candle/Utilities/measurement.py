from tqdm import tqdm

from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager

from slacker import sendslack
from Utilities.data import *

import seaborn as sns
   

def plot_IQ_generic(I, Q, n, xdata = [0.0],  title = "", xlabel = ""):
        
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    fig.suptitle(title, fontsize=11)
    
    ax1.plot(xdata, mag(I, Q))
    ax1.set_ylabel("mag (V)")
    
    ax2.plot(xdata, phase(I, Q))
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("phase (rad)")
    
    plt.tight_layout()
    
    return fig;


"""
getIQ(
      config,       :: configuration dictionary (usually imported from a file named config.py),
                       which is interpretted by QUA to determine settings on OPX 
      jobtype):     :: QUA jobtype (e.g. qubit_spec, rr_spec)
"""
def getIQ(config, jobtype,      showprogress = False, 
                                n_avg = 1000, 
                                notify = False, 
                                liveplot = False, 
                                plotfunc = plot_IQ_generic,
                                **plotkwargs):
    
    
    # Open Communication with the Server
    qmm = QuantumMachinesManager()

    # execute the job
    qm = qmm.open_qm(config)
    job = qm.execute(jobtype)
    res_handles = job.result_handles

    # get handles for I and Q
    I_handle = res_handles.get("I")
    Q_handle = res_handles.get("Q")
    n_handle = res_handles.get('n')
    I_handle.wait_for_values(1)
    Q_handle.wait_for_values(1)
    n_handle.wait_for_values(1)
    
    # display progress bar and send slack notification
    if showprogress:
        
        # get handle for iteration number n, and initialize counter
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
                    
    if liveplot:
        # make plot (initialize axes)
        
        while(I_handle.is_processing()):
            I = I_handle.fetch_all()
            Q = I_handle.fetch_all()
            n = n_handle.fetch_all()
            plotfunc(I, Q, n, **plotkwargs)
    else:
        res_handles.wait_for_all_values()
            
                    
    # retrieve I and Q values      
    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    
    # close quantum machine
    qmm.close_all_quantum_machines()
    
    return I, Q, job;


def plot_IQ_blobs_init():
    
    plot = sns.jointplot()
    plot.set_axis_labels('I', 'Q')
    plot.ax_marg_x.grid('on')
    plot.ax_marg_y.grid('on')
    plot.fig.tight_layout()
    ax = plt.gca()
    return plot, ax


def plot_IQ_blobs(plot, ax, datadict):
    
    I = datadict["I"]
    Q = datadict["Q"]
    I_exc = datadict["I_exc"]
    Q_exc = datadict["Q_exc"]
    
    data = {}
    data["I"] = I
    data["Q"] = Q
    
    plot.ax_joint.plot(I*1e3, Q*1e3, 'o', label = "ground")
    plot.ax_joint.plot(I_exc*1e3, Q_exc*1e3, 'o', label = "excited")
    # plot.plot_marginals(sns.kdeplot)
    ax.legend()
    


def get_results(config, jobtype, result_names = ["I", "Q", "n"],
                                showprogress = False, 
                                n_avg = 1000, 
                                notify = False, 
                                liveplot = False, 
                                plotfunc = plot_IQ_blobs,
                                plotinit = plot_IQ_blobs_init,
                                **plotkwargs):
    
    
    # Open Communication with the Server
    qmm = QuantumMachinesManager()

    # execute the job and get result handles
    qm = qmm.open_qm(config)
    job = qm.execute(jobtype)
    res_handles = job.result_handles
    
    
    # create a dictionary with result_names as keys
    # and their corresponding handles as values
    handles_dict = {}
    for name in result_names:
        handles_dict[name] = res_handles.get(name)
    
    # wait for first values to come in before continuing to run python code
    for handle in handles_dict.values():
        handle.wait_for_values(2)
        is_processing = lambda: handle.is_processing()
    
    # make progressmeter, if desired
    if showprogress:
        make_progress_meter(handles_dict["n"], n_avg)
    
    # make live plot, if desired
    if liveplot:
        # make plot (initialize axes)
        plot, ax = plotinit()
        
        while(is_processing()):
            datadict = get_data_from_handles(handles_dict)
            plotfunc(plot, ax, datadict, **plotkwargs)
    
    res_handles.wait_for_all_values()
            
                    
    # retrieve all values  
    datadict = get_data_from_handles(handles_dict)
    plotfunc(plot, ax, datadict)
    
    # close quantum machine
    qmm.close_all_quantum_machines()
    
    return datadict, job;





def get_data_from_handles(handles_dict):
    
    datadict = {}
    
    for (key, handle) in handles_dict.items():
        datadict[key] = handle.fetch_all()['value']
    
    return datadict




# display progress bar and send slack notification
def make_progress_meter(n_handle, n_total):
    
    # initialize counter
    n0 = 0 
    
    # create progressbar using tqdm
    with tqdm(total = n_total) as progress_bar:
        
        while(n_handle.is_processing()):
            
            
            n = n_handle.fetch_all()['value'][-1] # retrieve iteration value n
            Δn = n - n0 # calculate change in n since last update 

            if Δn > 0:
                progress_bar.update(Δn) # update progressbar with increase in n
                n0 = n # reset counter
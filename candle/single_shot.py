from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from config_MUN11 import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from file_utility import *
from plot_utiliity import *

def main(qm, dformat, bPlot = True, **kwargs):
    n_reps = 2000


    ###################
    # The QUA program #
    ###################
    
    with program() as IQ_blobs:
    
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        I_st_exc = declare_stream()
        Q_st_exc = declare_stream()
    
        resetTime=int(400e3)
        with for_(n, 0, n < n_reps, n + 1):
            
                wait(resetTime, "qubit")
                align("qubit", "rr") 
                measure("readout", "rr", None, *res_demod(I, Q))
                save(I, I_st)
                save(Q, Q_st)
                
                align('qubit','rr')
                
                wait(resetTime, "qubit")
                play("pi", "qubit")
                align("qubit", "rr")
                measure("readout", "rr", None, *res_demod(I, Q))
                save(I, I_st_exc)
                save(Q, Q_st_exc)
                
        with stream_processing():
            I_st.save_all('I')
            Q_st.save_all('Q')
            I_st_exc.save_all('I_exc')
            Q_st_exc.save_all('Q_exc')
        
    job = qm.execute(IQ_blobs)
    res_handles = job.result_handles
    res_handles.wait_for_all_values()
    I = res_handles.get("I").fetch_all()['value']
    Q = res_handles.get("Q").fetch_all()['value']
    I_exc = res_handles.get("I_exc").fetch_all()['value']
    Q_exc = res_handles.get("Q_exc").fetch_all()['value']
    
   
    # plt.figure()
    # plt.plot(I*1e3,Q*1e3,'.',alpha=0.5)
    # plt.axis('equal')
    # plt.plot(I_exc*1e3,Q_exc*1e3,'.',alpha=0.5)
    # plt.title('I_exc_mean=%.2f, '%(np.mean(I_exc*1e3))+'I_exc_std=%.2f,\\ '%(np.std(I_exc*1e3)) +
    #           'I_mean=%.2f, '%(np.mean(I*1e3)) +'I_std=%.2f, '%(np.std(I*1e3)) )
    
    dataDict = {'I': I,
                 'Q': Q,
                 'I_exc': I_exc,
                 'Q_exc': Q_exc}
    if dformat == 'csv':
        dataDf = pd.DataFrame(dataDict)
        datapath = filepath + datestr
        filename = 'IQblobs_wait=%dcycles_%s'%(resetTime,datestr)
        filename = make_new_filename(datapath, filename, fformat='csv')
        dataDf.to_csv(datapath + '\\' + filename)
        if bPlot:
            df, plot = plot_iq_hist(datapath + '\\' + filename)
            title = plot.fig._suptitle.get_text()
            title += '\n' + ' wait time = %.3f '%(resetTime * 4 /1e3)+ r'$\mu$' + 's'
            plot.fig.suptitle(title)
            plot.fig.savefig(datapath + '\\' + filename.replace('csv', 'png'))
            
    elif dformat == 'hdf5':
        if 'dname' in kwargs:
            dname = kwargs.get('dname')
        if 'dindex' in kwargs:
            dindex = kwargs.get('dindex')
        df = h5py.File(dname, 'a')
        if 'iq_blobs/%d' in df.keys():
            raise ValueError('iq_blobs/dindex = %d exists'%(dindex))
        else:
            gp = df.create_group('iq_blobs/%d'%(dindex))
            for key, value in dataDict.items():
                gp.create_dataset(key, data=value)
        
        df.close()

            
            

if __name__ == '__main__':
    ######################################
    # Open Communication with the Server #
    ######################################
    qmm = QuantumMachinesManager()
    
    ####################
    # Simulate Program #
    ####################
    # simulation_config = SimulationConfig(
    #                     duration=5000,
    #                     simulation_interface=LoopbackInterface([("con1", 3, "con1", 1), ("con1", 4, "con1", 2)]))
    
    # job = qmm.simulate(config, IQ_blobs, simulation_config)
    qm = qmm.open_qm(config)
    main(qm, dformat='csv', bPlot=True)
    




 
def plot_iq_hist(df):
    # [states.append(r'$|g\rangle$') for i in range(len(df['Q']))]
    # [states.append(r'$|e\rangle$') for i in range(len(df['Q']))]
    # i_threshold = np.mean(np.hstack((df['I']*1e3, df['I_exc']*1e3)))
    # print('i_threshold = %.2f mV'%(i_threshold))

    data = {
            'I (mV)':   np.hstack((df['I']*1e3, df['I_exc']*1e3)),
            'Q (mV)':   np.hstack((df['Q']*1e3, df['Q_exc']*1e3)),
            'states': states}
    dataF = pd.DataFrame(data=data,dtype=np.float32)
    plot = sns.jointplot(data=dataF, x='I (mV)',y='Q (mV)', hue='states')
    plot.ax_joint.axvline(x=i_threshold)
    plot.ax_marg_x.axvline(x=i_threshold)
    plot.fig.suptitle('readout fildeity = %.2f'%(1-p01-p10) + 
                      '\n'+'i threshold = %.2f mV'%(i_threshold))
    plot.fig.tight_layout()
    return df, plot



    
    


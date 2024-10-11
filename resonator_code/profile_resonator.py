# -*- coding: utf-8 -*-
"""
=======================================================
Program : resonator_profiling/profile_resonator.py
=======================================================
Summary:
"""
__author__ =  ["James Farmer", "Sadman Ahmed Shanto"]
__date__ = "10/28/2022"
__email__ = "shanto@usc.edu"

#libraries used
import Labber
import numpy as np
import os
import matplotlib.pyplot as plt
from fitTools.Resonator import Resonator
import logging

if __name__ == "__main__":
    plt.rcParams.update({'font.size':14})

    fpath = input("Path to .hdf5 file: ").replace('"','')
    power_offset = int(input("What is the total attenuation after the VNA? "))
    res_type = input("What type of resonator is this? (e.g. 'n','t','r') ")
    path,fname = os.path.split(fpath)
    path += r'\\'
    figpath = 'figures\\'+fname[:-4]+'\\'
    if not os.path.exists(path+'figures\\'):
        os.mkdir(path+'figures\\')
    if not os.path.exists(path+figpath):
        os.mkdir(path+figpath)
    lf = Labber.LogFile(path + fname)

    logFileName = path + f"profile_info_{fname[:-4]}log"
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler(logFileName, 'a'))
    print = logger.info

    nEntries = lf.getNumberOfEntries()
    try:
        power = np.squeeze(np.round(lf.getData(name='Agilent Network Analyzer - Output power'),decimals=2))
    except:
        power = -np.squeeze(np.round(lf.getData(name='Vaunix Lab Brick Digital Attenuator - Attenuation'),decimals=2))

    print(power)
    fits = {'f':[],'Q':[],'Qint':[],'Qext':[]}

    for n in range(nEntries):
        (xdata,ydata) = lf.getTraceXY(entry = n)
        res = Resonator(res_type,xdata,ydata)
        res.autofit()
        fits['f'].append(res.f0*1e-9)
        fits['Q'].append(res.Q)
        fits['Qint'].append(res.Qi)
        fits['Qext'].append(res.Qc)

        if res.fit_found:
            print(20*"=")
            # res.show(savefile = path+figpath+fname[:-4]+'resonance_{}-dBm.png'.format(power[n]))
            print('\nFit at {} dBm'.format(power[n]))
            print(res)
            print(20*"="+"\n")

    power += power_offset

    fig = plt.figure(figsize=[9,6],constrained_layout=True)
    plt.plot(power,fits['f'],'r.')
    plt.ylim([min(fits['f'])-0.1,max(fits['f'])+0.1])
    plt.title('Frequency vs power')
    plt.xlabel('Input Power [dBm]')
    #plt.xlabel('LO attenuation [dB]')
    plt.ylabel('Frequency from fit [GHz]')
    plt.savefig(path+figpath+fname[:-4]+'_f0-vs-P.png')
    plt.show()

    fig = plt.figure(figsize=[9,6],constrained_layout=True)
    plt.scatter(power,fits['Q'],s=20,c='r',label='Total Q')
    plt.scatter(power,fits['Qint'],s=20,c='b',label='Internal Q')
    plt.scatter(power,fits['Qext'],s=20,c='g',label='external Q')
    plt.yscale("log")
    plt.title('Q vs power')
    plt.xlabel('Input Power [dBm]')
    #plt.xlabel('LO attenuation [dB]')
    plt.ylabel('Quality factor')
    plt.legend()
    plt.savefig(path+figpath+fname[:-4]+'_Q-vs-P.png')
    plt.show()

    plt.scatter(power,np.array(fits['Q'])/np.array(fits['Qint']))
    plt.ylabel("Q/Q_i")
    plt.xlabel('Input Power [dBm]')
    plt.savefig(path+figpath+fname[:-4]+'_Q_over_Qi.png')
    plt.show()


    print(f"\n\nLog File: {logFileName}")
    print(f"\nPlots Directory: {path + figpath}")







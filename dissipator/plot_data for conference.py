# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 10:10:17 2023

@author: lfl
"""
import matplotlib.pyplot as plt
import numpy as np
import plot_functions as pf
import scienceplots
plt.style.use(["science",'no-latex'])
plt.rcParams["figure.dpi"] = 900
file1='G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\cavity-reset\\data_018.csv'
data1=np.genfromtxt(file1, delimiter=',', skip_header=3)
#print(data)
#np.shape(data)
t_arr=data1[0,:]
I1=data1[1,:]
Q1=data1[2,:]
ydata1 = np.abs(I1+1j*Q1)
fits, error = pf.fit_data(t_arr,ydata1,sequence='cavity-reset',dt=t_arr[-1]*1e-6/len(t_arr))
fitted_vals1=fits[0]*np.exp(-t_arr/fits[1])+fits[2]

file2='G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\cavity-reset\\data_019.csv'
data2=np.genfromtxt(file2, delimiter=',', skip_header=3)
#print(data)
#np.shape(data)
t_arr2=data2[0,:]
I2=data2[1,:]
Q2=data2[2,:]
ydata2 = np.abs(I2+1j*Q2)
fits2, error2 = pf.fit_data(t_arr2,ydata2,sequence='cavity-reset',dt=t_arr2[-1]*1e-6/len(t_arr2))
fitted_vals2=fits2[0]*np.exp(-t_arr2/fits2[1])+fits2[2]
fig, ax = plt.subplots(figsize=(6,6))
plt.scatter(t_arr,ydata1*1e3, label='with ffl, T2=%.2f $\pm$ %.2f $\mu$s'%(fits[1], error[1]), marker='x')
plt.plot(t_arr,fitted_vals1)
plt.scatter(t_arr2,ydata2*1e3, label='without ffl, T2=%.2f $\pm$ %.2f $\mu$s'%(fits2[1], error2[1]), marker='o')
plt.plot(t_arr2,fitted_vals2)
plt.xlabel("Time($\mu$s)")
plt.ylabel("Digitizer Voltage (mV)")
plt.title('Cavity Reset (wait time 500ns)')
ax.set_ylim([12.4, 13.4])
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig('G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\cavity-reset\\cavity reset.eps', format='eps')







dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\'
filename = 'cavity_ringdown_flux_ffl_len=0.0mA_fflFreq=3.05GHz_DA=30dB_fDA=13dB_rrLen=1000clks_navg=100000_4.h5'
f = h5py.File(dataDir + filename, 'r')
list(f[list(f.keys())[0]].keys())
for i, flux in enumerate(list(f.keys())):
    
    for j, amp in enumerate(list(f[flux].keys())):
        
        I=f[flux][amp]['I'][()]
        Q=f[flux][amp]['Q'][()]
        ydata = np.abs(I+1j*Q)
        t_arr=f[flux][amp]['time'][()]
        fits, error = pf.fit_data(t_arr,ydata,sequence='ringdown',dt=t_arr[-1]*1e-6/len(t_arr))
        fitted_vals3=fits[0]*np.exp(-t_arr/fits[1])+fits[2]
        #fig = pf.plot_data(t_arr,ydata,sequence='cavity-reset',fitted_pars=fitted_pars,nAverages=2000, pi2Width=qb.pars['pi_half_len'],
                         #qubitDriveFreq=qb.pars['qubit_LO']+qb.pars['qubit_IF'],qb_power = -8,iteration=0,amp_ffl_scale=j*0.1, amp=0.5,flux=0.7, error=error, ffl_len=ffl_len_list[i])

dataDir1 = 'G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230822\\'
filename1 = 'cavity_ringdown_flux_ffl_len=0.0mA_fflFreq=3.05GHz_DA=30dB_fDA=13dB_rrLen=1000clks_navg=50000_2.h5'
f1 = h5py.File(dataDir1 + filename1, 'r')
list(f1[list(f1.keys())[0]].keys())
for i, flux in enumerate(list(f1.keys())):
    
    for j, amp in enumerate(list(f1[flux].keys())):
        
        I=f1[flux][amp]['I'][()]
        Q=f1[flux][amp]['Q'][()]
        ydata1 = np.abs(I+1j*Q)
        t_arr1=f1[flux][amp]['time'][()]
        fits1, error1 = pf.fit_data(t_arr1,ydata1,sequence='ringdown',dt=t_arr1[-1]*1e-6/len(t_arr1))
        fitted_vals1=fits1[0]*np.exp(-t_arr1/fits1[1])+fits1[2]

fig, ax = plt.subplots(figsize=(6,6))
plt.scatter(t_arr,ydata*1e3, marker='x', label='without ffl, $T_{ring}$=%.3f $\pm$ %.4f $\mu$s'%(fits[1], error[1]))
plt.plot(t_arr,fitted_vals3)
plt.scatter(t_arr1,ydata1*1e3, label='with ffl, $T_{ring}$=%.2f $\pm$ %.3f $\mu$s'%(fits1[1], error1[1]), marker='o')
plt.plot(t_arr1,fitted_vals1)
plt.xlabel("Time($\mu$s)")
plt.ylabel("Digitizer Voltage (mV)")
plt.title('Cavity Ringdown')
plt.legend()
plt.tight_layout()
#ax.set_ylim([12.4, 13.4])
#plt.show()
plt.savefig('G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\cavity-reset\\cavity ringdown.eps', format='eps')
plt.show()



file1='G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230822\\cavity-cooling\\data_046.csv'
data1=np.genfromtxt(file1, delimiter=',', skip_header=3)
#print(data)
#np.shape(data)
t_arr=data1[0,:]
I1=data1[1,:]
Q1=data1[2,:]
ydata1 = np.abs(I1+1j*Q1)
fits, error = pf.fit_data(t_arr,ydata1,sequence='cavity-cooling',dt=t_arr[-1]*1e-6/len(t_arr))
fitted_vals1=fits[0]*np.exp(-t_arr/fits[1])+fits[2]

file2='G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230822\\cavity-cooling\\data_043.csv'
data2=np.genfromtxt(file2, delimiter=',', skip_header=3)
#print(data)
#np.shape(data)
t_arr2=data2[0,:]
I2=data2[1,:]
Q2=data2[2,:]
ydata2 = np.abs(I2+1j*Q2)
fits2, error2 = pf.fit_data(t_arr2,ydata2,sequence='cavity-cooling',dt=t_arr2[-1]*1e-6/len(t_arr2))
fitted_vals2=fits2[0]*np.exp(-t_arr2/fits2[1])+fits2[2]
fig, ax = plt.subplots(figsize=(8,6))
plt.scatter(t_arr,ydata1*1e3, label='with ffl, T2=%.2f $\pm$ %.2f $\mu$s'%(fits[1], error[1]), marker='x')
plt.plot(t_arr,fitted_vals1)
plt.scatter(t_arr2,ydata2*1e3, label='without ffl, T2=%.2f $\pm$ %.2f $\mu$s'%(fits2[1], error2[1]), marker='o', alpha=.3, )
plt.plot(t_arr2,fitted_vals2)
plt.xlabel("Time($\mu$s)")
plt.ylabel("Digitizer Voltage (mV)")
plt.title('Cavity Cooling (Improvement in bare T2)')
#ax.set_ylim([12.4, 13.4])
plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig('G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\cavity-reset\\cavity cooling1.eps', format='eps')


file1='G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230823\\cavity-cooling\\data_025.csv'
data1=np.genfromtxt(file1, delimiter=',', skip_header=3)
#print(data)
#np.shape(data)
t_arr=data1[0,:]
I1=data1[1,:]
Q1=data1[2,:]
ydata1 = np.abs(I1+1j*Q1)
fits, error = pf.fit_data(t_arr,ydata1,sequence='cavity-cooling',dt=t_arr[-1]*1e-6/len(t_arr))
fitted_vals1=fits[0]*np.exp(-t_arr/fits[1])+fits[2]

file2='G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230823\\cavity-cooling\\data_034.csv'
data2=np.genfromtxt(file2, delimiter=',', skip_header=3)
#print(data)
#np.shape(data)
t_arr2=data2[0,:]
I2=data2[1,:]
Q2=data2[2,:]
ydata2 = np.abs(I2+1j*Q2)
fits2, error2 = pf.fit_data(t_arr2,ydata2,sequence='cavity-cooling',dt=t_arr2[-1]*1e-6/len(t_arr2))
fitted_vals2=fits2[0]*np.exp(-t_arr2/fits2[1])+fits2[2]
fig, ax = plt.subplots(figsize=(8,6))
plt.scatter(t_arr,ydata1*1e3, label='without ffl, T2=%.2f $\pm$ %.2f $\mu$s'%(fits[1], error[1]), marker='x')
plt.plot(t_arr,fitted_vals1)
plt.scatter(t_arr2,ydata2*1e3, label='with ffl', marker='o', alpha=.3, )
plt.plot(t_arr2,fitted_vals2)
plt.xlabel("Time($\mu$s)")
plt.ylabel("Digitizer Voltage (mV)")
plt.title('Cavity Cooling (with cavity drive)')
#ax.set_ylim([12.4, 13.4])
plt.legend()
plt.tight_layout()
plt.show()
plt.savefig('G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\cavity-reset\\cavity cooling with cavity.eps', format='eps')



data=np.load('G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\2d sweep cavity cooling new.npy')
fig, ax = plt.subplots(figsize=(6,6))
#plt.rcParams["figure.dpi"] = 900
plot=ax.imshow(data, interpolation='nearest',  extent=[0,0.35,0,0.26], vmin=0.3, vmax=8, aspect="auto", origin="lower")
#ax.set_yticklabels([500,1000,1800,2800,4000])
plt.xlabel("ffl amp")
plt.ylabel("cavity drive amp")
plt.title('Cavity Cooling')
plt.tight_layout()
plt.colorbar(plot)
plt.savefig('G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230821\\cavity-reset\\cavity cooling 2d sweep.eps', format='eps')

for i in [5,6,7,8,9,10]:
    plt.plot(np.linspace(0,0.35,11)**2,1/(data[:,i]))
plt.xlabel("cavity drive power")
plt.ylabel("dephaing rate (1/us)")
plt.tight_layout()
plt.show()




dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230913\\interleaved_exp\\'
filename = 'interleaved_cavity_reset_sweepffllen_flux=1150.0mA_fflFreq=2.98GHz_DA=29dB_fDA=9dB_amp_ffl_scale=1.0_navg=15000amp_r_scale=0.8_2.h5'
f = h5py.File(dataDir + filename, 'r')
list(f[list(f.keys())[0]].keys())

f.keys()
tb=np.zeros(20)
tc=np.zeros(20)
tcf=np.zeros(20)
eb=np.zeros(20)
ec=np.zeros(20)
ecf=np.zeros(20)


for j in list(f.keys()):
    for i, ffl_len in enumerate(list(f[j].keys())):
        yb=np.abs(f[j][ffl_len]['I_b'][()]+1j*f[j][ffl_len]['Q_b'][()])
        yc=np.abs(f[j][ffl_len]['I_c'][()]+1j*f[j][ffl_len]['Q_c'][()])
        ycf=np.abs(f[j][ffl_len]['I_cf'][()]+1j*f[j][ffl_len]['Q_cf'][()])
        tarr=f[j][ffl_len]['t_arr'][()]
        fitsb, errorb = pf.fit_data(tarr,yb,sequence='cavity-reset',dt=tarr[-1]*1e-6/len(tarr))
        fitsc, errorc = pf.fit_data(tarr,yc,sequence='cavity-reset',dt=tarr[-1]*1e-6/len(tarr))
        fitscf, errorcf = pf.fit_data(tarr,ycf,sequence='cavity-reset',dt=tarr[-1]*1e-6/len(tarr))
        #print(fitsb)
        #tb[i]=fitsb[1]
        #tc[i]=fitsc[1]
        #tcf[i]=fitscf[1]
        #tc[tc>7]=0.1
        #tc[tc<0]=0.1
        eb[i]=errorb[1]
        ec[i]=errorc[1]
        ecf[i]=errorcf[1]
        tb[i]=f[j][ffl_len]['fitted_pars_b'][()][1]
        tc[i]=f[j][ffl_len]['fitted_pars_c'][()][1]
        tcf[i]=f[j][ffl_len]['fitted_pars_cf'][()][1]

ec[ec>0.2]=0.
t_arr=np.array([1004,128,1312,1524,164,1736,1948,196,2156,232,2368,2580,260,2788,3000,472,60,680,892,96]) +100     
fig, ax = plt.subplots(figsize=(8,6))
plt.errorbar(t_arr,tb, fmt='x', label='bare T2', yerr=eb)
plt.errorbar(t_arr,tc, fmt='o', label='T2 with cavity', yerr=ec)
plt.errorbar(t_arr,tcf,  fmt='+', label= 'T2 with ffl and cavity', yerr=ecf)
plt.xlabel("ffl length (ns)")
plt.ylabel("T2 (in us)")
plt.title('Interleaved cavity reset(amp 0.8, ffl atten 9)')
plt.legend()
plt.tight_layout()




dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230914\\interleaved_exp\\'
filename = 'interleaved_cavity_cooling_sweepffllen_flux=1150.0mA_fflFreq=2.98GHz_DA=29dB_fDA=9dB_amp_ffl_scale=1.0_navg=15000amp_r_scale=0.15_3.h5'
f = h5py.File(dataDir + filename, 'r')
list(f[list(f.keys())[0]].keys())

f.keys()
tb=np.zeros(7)
tc=np.zeros(7)
tcf=np.zeros(7)
eb=np.zeros(7)
ec=np.zeros(7)
ecf=np.zeros(7)


for j in list(f.keys()):
    for i, atten in enumerate(list(f[j].keys())):
        print(atten)
        yb=np.abs(f[j][atten]['I_b'][()]+1j*f[j][atten]['Q_b'][()])
        yc=np.abs(f[j][atten]['I_c'][()]+1j*f[j][atten]['Q_c'][()])
        ycf=np.abs(f[j][atten]['I_cf'][()]+1j*f[j][atten]['Q_cf'][()])
        tarr=f[j][atten]['t_arr'][()]
        fitsb, errorb = pf.fit_data(tarr,yb,sequence='cavity-cooling',dt=tarr[-1]*1e-6/len(tarr))
        fitsc, errorc = pf.fit_data(tarr,yc,sequence='cavity-cooling',dt=tarr[-1]*1e-6/len(tarr))
        fitscf, errorcf = pf.fit_data(tarr,ycf,sequence='cavity-cooling',dt=tarr[-1]*1e-6/len(tarr))
        #print(fitsb)
        #tb[i]=fitsb[1]
        #tc[i]=fitsc[1]
        #tcf[i]=fitscf[1]
        #tc[tc>7]=0.1
        #tc[tc<0]=0.1
        eb[i]=errorb[1]
        ec[i]=errorc[1]
        ecf[i]=errorcf[1]
        tb[i]=f[j][atten]['fitted_pars_b'][()][1]
        tc[i]=f[j][atten]['fitted_pars_c'][()][1]
        tcf[i]=f[j][atten]['fitted_pars_cf'][()][1]

ec[ec>0.2]=0.
ecf[ecf>0.2]=0.
atten_list=np.array([11,13,15,17,20,23,9])  
fig, ax = plt.subplots(figsize=(8,6))
plt.errorbar(atten_list,tb, fmt='x', label='bare T2', yerr=eb)
plt.errorbar(atten_list,tc, fmt='o', label='T2 with cavity', yerr=ec)
plt.errorbar(atten_list,tcf,  fmt='+', label= 'T2 with ffl and cavity', yerr=ecf)
ax.invert_xaxis()
plt.xlabel("atten (db)")
plt.ylabel("T2 (in us)")
plt.title('Interleaved cavity cooling(cavity amp 0.15)')
plt.legend()
plt.tight_layout()


T2=np.zeros((11,11))
err=np.zeros((11,11))
dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230916\\interleaved_exp\\'
filename = 'interleaved_cavity_cooling_sweepffllen_flux=1150.0mA_fflFreq=2.98GHz_DA=29dB_fDA=9dB_amp_ffl_scale=1.0_navg=12000amp_r_scale=1_1.h5'
f = h5py.File(dataDir + filename, 'r')
list(f[list(f.keys())[0]].keys())
for j in list(f.keys()):
    for i, amp in enumerate(list(f[j].keys())):
        for k, atten in enumerate(list(f[j][amp])):
            T2[i,k]=f[j][amp][atten]['fitted_pars'][()][1]
            err[i,k]=f[j][amp][atten]['error'][()][1]
T2[T2<0]=1.6
err[err>1]=0.
fig, ax = plt.subplots(figsize=(6,6))
#plt.rcParams["figure.dpi"] = 900
plot=ax.imshow(T2, interpolation='nearest',  extent=[9,30,0,0.2], vmin=0.3, vmax=8, aspect="auto", origin="lower")
#ax.set_yticklabels([500,1000,1800,2800,4000])
plt.xlabel("ffl atten")
plt.ylabel("cavity drive amp")
plt.title('Cavity Cooling')
plt.tight_layout()
plt.colorbar(plot)


fig, ax = plt.subplots(figsize=(8,6))
for i, atten in enumerate([9, 10, 11, 12, 13, 15, 17, 19, 23, 26, 30]):
    plt.errorbar(np.linspace(0,0.2,11)**2,1/T2[:,i], label='fixed ffl atten = %d'%(atten), yerr=np.divide(err[:,i], (T2[:,i]**2)))
plt.xlabel("cav power ($scale^2$)")
plt.ylabel("dephasing rate (in $us^{-1}$)")
#ax.invert_xaxis()
plt.title('cavity cooling')
plt.legend()
plt.tight_layout()


T2low=np.zeros((11,6))
errlow=np.zeros((11,6))
dataDir = 'G:\\Shared drives\\CavityCooling\\data\\diss08_11a\\20230917\\interleaved_exp\\'
filename = 'interleaved_cavity_cooling_sweepffllen_flux=1150.0mA_fflFreq=2.98GHz_DA=29dB_fDA=9dB_amp_ffl_scale=1.0_navg=12000amp_r_scale=1_2.h5'
f = h5py.File(dataDir + filename, 'r')
list(f[list(f.keys())[0]].keys())
for j in list(f.keys()):
    for i, amp in enumerate(list(f[j].keys())):
        for k, atten in enumerate(list(f[j][amp])):
            T2low[i,k]=f[j][amp][atten]['fitted_pars'][()][1]
            errlow[i,k]=f[j][amp][atten]['error'][()][1]
T2[T2<0]=1.6
err[err>1]=0.
fig, ax = plt.subplots(figsize=(6,6))
#plt.rcParams["figure.dpi"] = 900
plot=ax.imshow(T2low, interpolation='nearest',  extent=[4,9,0,0.2], vmin=0.3, vmax=8, aspect="auto", origin="lower")
#ax.set_yticklabels([500,1000,1800,2800,4000])
plt.xlabel("ffl atten")
plt.ylabel("cavity drive amp")
plt.title('Cavity Cooling')
plt.tight_layout()
plt.colorbar(plot)
pf.heatplot(xdata=[9,10, 11, 12, 13, 15, 17, 19, 23, 26, 30], ydata=np.linspace(0,0.2,11), data=T2, xlabel="ffl atten", ylabel="cavity drive amp", cbar_label='T2 (us)')

T2new=np.delete(T2,0,1)
t2f=np.append(T2low,T2new, axis=1)
fig, ax = plt.subplots(figsize=(6,6))
#plt.rcParams["figure.dpi"] = 900
plot=ax.imshow(t2f, interpolation='none',  extent=[4,30,0,0.2], vmin=0.3, vmax=8, aspect="auto", origin="lower")
plt.xticks([4,5,6,7,8,9,10, 11, 12, 13, 15, 17, 19, 23, 26, 30])
#ax.set_yticklabels([500,1000,1800,2800,4000])
plt.xlabel("ffl atten")
plt.ylabel("cavity drive amp")
plt.title('Cavity Cooling')
plt.tight_layout()
plt.colorbar(plot)
pf.heatplot(xdata=[4,5,6,7,8,9,10, 11, 12, 13, 15, 17, 19, 23, 26, 30], ydata=np.linspace(0,0.2,11), data=t2f, xlabel="ffl atten", ylabel="cavity drive amp", cbar_label='T2 (us)')

fig, ax = plt.subplots(figsize=(8,6))
for i, atten in enumerate([9, 10, 11, 12, 13, 15, 17, 19, 23, 26, 30]):
    plt.errorbar(np.linspace(0,0.2,11)**2,1/T2[:,i], label='fixed ffl atten = %d'%(atten), yerr=np.divide(err[:,i], (T2[:,i]**2)))
plt.xlabel("cav power ($scale^2$)")
plt.ylabel("dephasing rate (in $us^{-1}$)")
#ax.invert_xaxis()
plt.title('cavity cooling')
plt.legend()
plt.tight_layout()
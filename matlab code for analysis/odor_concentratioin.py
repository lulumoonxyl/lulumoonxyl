# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 11:41:03 2022

@author: feihu
"""
#%% load library
import numpy as np
import pyabf
import matplotlib.pyplot as plt
from itertools import product
#%% get the data from trace
V=np.full((960000,3),np.nan)
t=np.full((960000,3),np.nan)
abf=pyabf.ABF('D:/Low suction, low gain 3_15min/H2O_row1_col2_0000.abf')
V[:,0]=abf.sweepY
t[:,0]=abf.sweepX

abf=pyabf.ABF('D:/Low suction, low gain 3_15min/H2O_row1_col2_0001.abf')
V[:,1]=abf.sweepY
t[:,1]=abf.sweepX

abf=pyabf.ABF('D:/Low suction, low gain 3_15min/H2O_row1_col2_0002.abf')
V[:,2]=abf.sweepY
t[:,2]=abf.sweepX
#%% plot the whole trace
Vmean1=np.mean(V[:2500],axis=0)
V1=V-Vmean1
plt.plot(t,V1,linewidth=0.5)
plt.xlabel('Time(s)')
plt.ylabel('Odor concentration (mV)')
plt.legend(['1st test', '2nd','3rd'])
plt.title('H2O_row1_col2')
# plt.savefig('C:/Users/feihu/OneDrive - McGill University/matlab/Analysis based on larvae/olfactory_output_pdf/odor gradient measurement/H2O_row1_col2_0000.svg')
#%%
## the recording is 16min long-->contain 960000 time frame
## record 6 times 

#10n2GA_row1_col2_0000,0001,0002
# t1=np.array([[202500,380000,562000,779000,923000], [203000,387000,552000,737000,922000],[209000,381000,560000,748000,915000]])
#10n2GA_row3_col2_0000,0001,0002
# t1=np.array([[206000,384000,565000,742500,927500], [199500,387500,570000,741000,930000],[193000,374500,544500,726000,907000]])
#H2O_row3_col2_0000,0001,0002
# t1=np.array([[165000,350000,535000,710000,884000], [129500,372500,548500,733000,910000],[131500,366500,557000,762500,922000]])
#H2O_row1_col2_0000,0001,0002
t1=np.array([[181000,370000,536800,725000,898500], [189200,372000,545000,731000,921500],[200000,371000,561000,737000,917000]])

#%%
Vlen=40000
V0=np.full((Vlen,6,3),np.nan)
t0=np.full((Vlen,6,3),np.nan)
Vmean=np.zeros((6,3))

for i, j in product(range(6),range(3)):

    if i==0:
        V0[:,i,j]=V[0:Vlen,j]
        t0[:,i,j]=t[0:Vlen,j]
        Vmean[i,j]=np.mean(V0[:2500,i,j])
    else: 
        len1=len(V[t1[j,i-1]:t1[j,i-1]+Vlen,j])
        V0[0:len1,i,j]=V[t1[j,i-1]:t1[j,i-1]+Vlen,j]
        t0[0:len1,i,j]=t[t1[j,i-1]:t1[j,i-1]+Vlen,j]
        Vmean[i,j]=np.mean(V0[:2500,i,j])

t0=t0-t0[0,:,:]
V0=V0-Vmean

#%% plt.subplot(1,3,1)
for i, j in product(range(6),range(3)):
    plt.subplot(3,1,j+1)
    plt.plot(t0[:,i,j],V0[:,i,j],linewidth=0.4)
    plt.ylabel('Odor [] (mV)')
    if j==2:
        plt.xlabel('Time(s)')
    elif j==0:
        plt.legend(['0min',' 3min','6min','9min','12min','15min'],fontsize=6, mode = "expand", ncol = 6)
        plt.title('H2O_row1_col2')
    plt.ylim(-1, 7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.grid(visible=True)
plt.savefig('C:/Users/feihu/OneDrive - McGill University/matlab/Analysis based on larvae/olfactory_output_pdf/odor gradient measurement/H2O_row1_col2.svg')


 
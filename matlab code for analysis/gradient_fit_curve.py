# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 14:59:23 2022

@author: feihu
"""
#%% import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#%% load data
data=pd.read_excel('C:/Users/Fei Huang/OneDrive - McGill University/matlab/gradient_areana_measurement.xlsx',
                 sheet_name='Total Plot')
df=pd.DataFrame(data,columns=['C13','C14','C15','C16','C17'],dtype=float)
gradient=df.iloc[0:11]
#%% fit a POLY curve onto different conditions (C13-C16)
colors=['orange','purple','yellow','blue','green']
#pos=np.arange(1,12)
pos=np.linspace(7.5,217.5,11) # x data
polyline = np.linspace(0,250,960) # x for fitted curve later
all_model=[]
all_coeff=np.zeros((5,4))
#%% find the parameter value for the fitted curve
def sigmoid(x,L,x0,k,b):
    y=L/(1+np.exp(-k*(x-x0)))+b
    return (y)

for i in range(5):
    if i!=4: 
        model=np.poly1d(np.polyfit(pos,gradient.iloc[:,i],3))
        plt.plot(polyline,model(polyline),color=colors[i])       
        plt.scatter(pos,gradient.iloc[:,i],color=colors[i])
        all_model.append(model)
        all_coeff[i,:]=model.coeffs
        print(model)
    else: 
        df_C17=gradient.iloc[:,i] # fit a sigmoid for C17
        p0=[max(df_C17),np.median(pos),1,min(df_C17)]
        popt,pcov=curve_fit(sigmoid,pos,df_C17,p0,method='dogbox')
        y=sigmoid(polyline,*popt)
        plt.plot(polyline,y,color=colors[i])
        plt.scatter(pos,gradient.iloc[:,i],color=colors[i])
        all_coeff[i,:]=popt

plt.xlabel('Position')
x=np.linspace(7.5,217.5,11,dtype=float)
#plt.xticks(pos,x)
plt.ylabel('Light Intensity (uW/mm^2)')
plt.legend(['C13','C14','C15','C16','C17'])
plt.title('The light intensity vs postions for different conditions')
plt.xlim([-10 ,250])
plt.savefig('C:/Users/feihu/OneDrive - McGill University/matlab/Analysis based on larvae/code/gradient turn analysis/light_gradient_avg.png')
plt.show()



#%% Use every single measurement instead of average
ind=pd.DataFrame(data,columns=['Position'])
ind=ind-1
df1=pd.DataFrame(data,columns=['C13.1','C14.1','C15.1','C16.1','C17.1'],dtype=float)
pos1=np.squeeze(x[ind])
#%%
all_mod=[]
polyline = np.linspace(0,250,960)
all_coef=np.zeros((5,4))
for i in range(5):
    if i!=4:
        model=np.poly1d(np.polyfit(pos1,df1.iloc[:,i],3))
        plt.plot(polyline,model(polyline),color=colors[i])       
        plt.scatter(pos1,df1.iloc[:,i],color=colors[i])
        all_mod.append(model)
        all_coef[i,:]=model.coeffs
        print(model)
    else:
        df_C17=df1.iloc[:,i] # fit a sigmoid for C17
        p0=[max(df_C17),np.median(pos),1,min(df_C17)]
        popt17,pcov17=curve_fit(sigmoid,pos1,df_C17,p0,method='dogbox')
        y=sigmoid(polyline,*popt17)
        plt.plot(polyline,y,color=colors[i])
        plt.scatter(pos1,df1.iloc[:,i],color=colors[i])
        all_coef[i,:]=popt17

plt.xlabel('Position')
#plt.xticks(pos,x)
plt.ylabel('Light Intensity (uW/mm^2)')
plt.legend(['C13','C14','C15','C16','C17'])
plt.title('The light intensity vs postions for different conditions')
plt.xlim([-10 ,250])
plt.savefig('C:/Users/feihu/OneDrive - McGill University/matlab/Analysis based on larvae/code/gradient turn analysis/light_gradient.png')
plt.show()
#%%
np.savetxt('C:/Users/feihu/OneDrive - McGill University/matlab/Analysis based on larvae/code/gradient turn analysis/light_gradient_coef.txt',all_coef,delimiter=' ')

#%% plot heatmap for each condition
import seaborn as sns
sns.heatmap(gradient.iloc[0:1,:],cmap='Oranges')
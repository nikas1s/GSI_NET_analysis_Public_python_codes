#This is a script to access the databses_complete and retrieve mass information for AME16 and AME20 and plots
#in a 2d colorplot with X=N and Y=Z space.

import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
###########################################################################################
#############           Subroutines                              ##########################
###########################################################################################

#This subroutine takes the chopped part of the dataframe and converts it in a flat 2d array
def convert_to_numpy_flat_array(array):
    array=array.to_numpy()
    array=array.flatten()
    return(array)

#This subroutine returns
def return_A_Z_mass_excess(df1,mass_model): #df1 is the complete database
    extrapolated=[]
    if mass_model=="mass_excess_AME16":
        extrapolated=df1[['extrapolated_AME16']][(df1[mass_model]!='nan')]
        error='uncertainty_AME16'
    if mass_model=="mass_excess_AME20":
        extrapolated=df1[['extrapolated_AME20']][(df1[mass_model]!='nan')]
        error='uncertainty_AME20'

    A=df1[['A']][(df1[mass_model]!="nan")]
    Z=df1[['Z']][(df1[mass_model]!="nan")]
    mass_excess=df1[[mass_model]][(df1[mass_model]!="nan")]
    unc=df1[[error]][(df1[error]!="nan")]

    A=convert_to_numpy_flat_array(A)
    Z=convert_to_numpy_flat_array(Z)
    extrapolated=convert_to_numpy_flat_array(extrapolated)
    mass_excess=convert_to_numpy_flat_array(mass_excess)
    unc=convert_to_numpy_flat_array(unc)
    #print (len(A),len(Z),len(mass_excess))
    #for i in range (0,len(Z)):
    #    print(A[i],Z[i],mass_excess[i])
    #print(extrapolated)
    return(A,Z,mass_excess,extrapolated,unc)


def simple_plot_initial_definitions():

    fig=plt.figure(figsize=(15,10))
    fig.set_size_inches(10, 6.5)
    font = {'family' : 'normal','size'   : 14.5}
    plt.rc('font', **font)
    return(ax)

def plot_stable_nuclei_as_black_squares():
    stable=np.loadtxt("/Users/stynikas/Data/Gabriel/Stable_Nuclides.txt",usecols=(0,1),unpack=True)
    St=np.empty((120,350))
    St[:]=np.nan
    for i in range (0, len(stable[1])):
        St[int(stable[1][i])][int(stable[0][i]-stable[1][i])]=1
    #plot the stable nuclei as black squares
    stable=ax.imshow(St,cmap="Greys_r")
    return(stable)


def plot_vlines(ax):
    print (i)
    #this is a test function
    return()

def box(Z,N):
    ax.hlines(y=Z-0.5,xmin=N-0.5,xmax=N+0.5,color='k')
    ax.vlines(x=N-0.5,ymin=Z-0.5,ymax=Z+0.5,color='k')
    ax.hlines(y=Z+0.5,xmin=N-0.5,xmax=N+0.5,color='k')
    ax.vlines(x=N+0.5,ymin=Z-0.5,ymax=Z+0.5,color='k')
    return()


def magic():
    magic=[20,28,50,82,126]
    limits_v_up=[20,25,40,62,100]
    limits_v_down=[8,14,25,45,70]
    limits_h_up=[40,50,90,126,126]
    limits_h_down=[20,30,60,126,126]
    j=0
    for i in magic:
        #print(limits_v_down[j])
        ax.vlines(x=i+0.5,ymin=limits_v_down[j],ymax=limits_v_up[j])
        ax.vlines(x=i-0.5,ymin=limits_v_down[j],ymax=limits_v_up[j])
        ax.hlines(y=i+0.5,xmin=limits_h_down[j],xmax=limits_h_up[j])
        ax.hlines(y=i-0.5,xmin=limits_h_down[j],xmax=limits_h_up[j])
        j=j+1
    return()

###########################################################################################
#############               Main Code                            ##########################
###########################################################################################

df1 = pd.read_pickle('/Users/stynikas/Python_codes/databases_codes/Database_complete_updated.pkl')


#############################################
#Mass excess nz plane plots HNPS2021
#############################################
#font = {'family' : 'normal','size'   : 14.5}
#plt.rc('font', **font)
n_z_plane=np.zeros((120,350))
n_z_plane[:]=np.nan
fig,ax=plt.subplots()
fig.set_size_inches(10, 4.5)

mass_model="mass_excess_AME16" #or mass_excess_FRDM mass_excess_AME20
A0,Z0,mass_excess16,extrapolated16,error16=return_A_Z_mass_excess(df1,mass_model)
mass_model="mass_excess_AME20" #or mass_excess_FRDM mass_excess_AME20
A1,Z1,mass_excess20,extrapolated20,error20=return_A_Z_mass_excess(df1,mass_model)
f = open("mass_excess",'w')
for i in range (0,len(A1)):
    f.write("{} {} {} {} \n".format(A1[i],Z1[i],mass_excess20[i],extrapolated20[i]))
Mass_excess_2d_array_diff=np.zeros((120,350))
Mass_excess_2d_array_diff[:]=np.nan
for i in range (0,len(A0)):
    for j in range (0,len(A1)):
        if A0[i]==A1[j] and Z0[i]== Z1[j]:
            if extrapolated16[i]=='no'and extrapolated20[j]=='no':
                Mass_excess_2d_array_diff[int(Z0[i])][int(A0[i]-Z0[i])]=mass_excess16[i]
                Mass_excess_2d_array_diff[int(Z1[j])][int(A1[j]-Z1[j])]=Mass_excess_2d_array_diff[int(Z1[j])][int(A1[j]-Z1[j])]-mass_excess20[j]

im=ax.imshow(Mass_excess_2d_array_diff,cmap="seismic",vmin=-0.5, vmax=0.5)
cbar=plt.colorbar(im)
stable=np.loadtxt("/Users/stynikas/Data/Gabriel/Stable_Nuclides.txt",usecols=(0,1),unpack=True)
St=np.empty((120,350))
St[:]=np.nan
for i in range (0, len(stable[1])):
    St[int(stable[1][i])][int(stable[0][i]-stable[1][i])]=1
stable=ax.imshow(St,cmap="Greys_r")
ax.set_ylim(8,80)
for i in range (0,len(A1)):
    if extrapolated20[i]=='yes':
        ax.plot(float(A1[i]-Z1[i]), Z1[i],marker='|',color='k')
for i in range (0,len(A0)):
    if extrapolated16[i]=='yes':
        ax.plot(float(A0[i]-Z0[i]), Z0[i],marker='_',color='k')
for i in range (0,len(A1)):
    if extrapolated20[i]=='no':
        ax.plot(float(A1[i]-Z1[i]), Z1[i],marker='x',color='k')
ax.set_xlim(10,120)
ax.set_ylabel('Z')
ax.set_xlabel('N')
magic=[20,28,50,82,126]
limits_v_up=[30,35,55,82,100]
limits_v_down=[8,10,20,40,70]
limits_h_up=[40,50,90,126,126]
limits_h_down=[10,20,50,126,126]
j=0
for i in magic:
    #print(limits_v_down[j])
    ax.vlines(x=i+0.5,ymin=limits_v_down[j],ymax=limits_v_up[j])
    ax.vlines(x=i-0.5,ymin=limits_v_down[j],ymax=limits_v_up[j])
    ax.hlines(y=i+0.5,xmin=limits_h_down[j],xmax=limits_h_up[j])
    ax.hlines(y=i-0.5,xmin=limits_h_down[j],xmax=limits_h_up[j])
    j=j+1
cbar.set_label(r'$M_{excess(AME16)}$-$M_{excess(AME20)}$ (MeV)')

plt.show()




#############################################
#Mass error nz plane plot
#############################################
fig,ax=plt.subplots()
fig.set_size_inches(10, 4.5)
mass_model="mass_excess_AME16" #or mass_excess_FRDM mass_excess_AME20
A0,Z0,mass_excess16,extrapolated16,error16=return_A_Z_mass_excess(df1,mass_model)
mass_model="mass_excess_AME20" #or mass_excess_FRDM mass_excess_AME20
A1,Z1,mass_excess20,extrapolated20,error20=return_A_Z_mass_excess(df1,mass_model)
f = open("mass_excess",'w')
for i in range (0,len(A1)):
    f.write("{} {} {} {} \n".format(A1[i],Z1[i],mass_excess20[i],extrapolated20[i]))
Mass_excess_2d_array_diff=np.zeros((120,350))
Mass_excess_2d_array_diff[:]=np.nan
betas=np.loadtxt('/Users/stynikas/Python_codes/beta-exp.d',unpack=True,skiprows=1)
for i in range (0,len(A0)):
    for j in range (0,len(A1)):
        if A0[i]==A1[j] and Z0[i]== Z1[j]:
            if extrapolated20[j]=='no':
                Mass_excess_2d_array_diff[int(Z0[i])][int(A0[i]-Z0[i])]=float(error20[j])
                if float(error20[j]>0.010):
                    for s in range (0,len(betas[1])):
                        if int(betas[0][s])==Z1[j] and int(betas[1][s])==(A1[j]-Z1[j]) and float(betas[2][s])>=0.10:
                            print (betas[:,s])
                            ax.plot(float(A1[j]-Z1[j]), Z1[j],marker='o',color='k',markersize=9,fillstyle='none')
im=ax.imshow(Mass_excess_2d_array_diff,cmap="Greens",vmin=0, vmax=0.1)
cbar=plt.colorbar(im)
stable=np.loadtxt("/Users/stynikas/Data/Gabriel/Stable_Nuclides.txt",usecols=(0,1),unpack=True)
St=np.empty((120,350))
St[:]=np.nan
for i in range (0, len(stable[1])):
    St[int(stable[1][i])][int(stable[0][i]-stable[1][i])]=1
    #plot the stable nuclei as black squares
stable=ax.imshow(St,cmap="Greys_r")
ax.set_ylim(8,80)
for i in range (0,len(A1)):
    if extrapolated20[i]=='yes':
        ax.plot(float(A1[i]-Z1[i]), Z1[i],marker='|',color='k')
for i in range (0,len(A0)):
    if extrapolated16[i]=='yes':
        ax.plot(float(A0[i]-Z0[i]), Z0[i],marker='_',color='k')
for i in range (0,len(A1)):
    if extrapolated20[i]=='no':
        ax.plot(float(A1[i]-Z1[i]), Z1[i],marker='x',color='k')
ax.set_xlim(10,120)
ax.set_ylabel('Z')
ax.set_xlabel('N')
magic=[20,28,50,82,126]
limits_v_up=[30,35,55,82,100]
limits_v_down=[8,10,20,40,70]
limits_h_up=[40,50,90,126,126]
limits_h_down=[10,20,50,126,126]
j=0
for i in magic:
    #print(limits_v_down[j])
    ax.vlines(x=i+0.5,ymin=limits_v_down[j],ymax=limits_v_up[j])
    ax.vlines(x=i-0.5,ymin=limits_v_down[j],ymax=limits_v_up[j])
    ax.hlines(y=i+0.5,xmin=limits_h_down[j],xmax=limits_h_up[j])
    ax.hlines(y=i-0.5,xmin=limits_h_down[j],xmax=limits_h_up[j])
    j=j+1
cbar.set_label(r'$M_{err (AME20)}$ (MeV)')
plt.show()

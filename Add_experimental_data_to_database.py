import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import fortranformat as ff
from os import listdir
import glob
from tabulate import tabulate
from scipy.optimize import curve_fit
import os
import math

#subroutine to calculate Q-values
def calculate_Qvalues(reaction_proj,reaction_product,name,i,j,df):
    #constants
    alpha_mass=4.00260325413
    neutron_mass=1.00866491588
    proton_mass=1.00782503224
    MeV=931.49432
    models=['I261','I261plusAME20']
    for model in models:
        name='{}{}_{}'.format(reaction_proj,reaction_product,model)
        mass='mass_{}'.format(model)
        reaction_unc='reaction_unc_{}{}_{}'.format(reaction_proj,reaction_product,model)
        try:
            if pd.isna(df[mass][i])==True:
                mass='mass_AME20'

            if reaction_proj == 'a':
                parent = (df[mass][i]+alpha_mass)*MeV
            elif reaction_proj=='n':
                parent = (df[mass][i]+neutron_mass)*MeV
            elif reaction_proj=='p':
                parent = (df[mass][i]+proton_mass)*MeV
            elif reaction_proj=='g':
                parent = (df[mass][i])*MeV

            mass='mass_{}'.format(model)
            if pd.isna(df[mass][j])==True:
                mass='mass_AME20'
            if reaction_product == 'a':
                daughter = (df[mass][j]+alpha_mass)*MeV
            elif reaction_product == 'n':
                daughter = (df[mass][j]+neutron_mass)*MeV
            elif reaction_product == '2n':
                daughter = (df[mass][j]+2.*neutron_mass)*MeV
            elif reaction_product == 'g':
                daughter = (df[mass][j])*MeV
            elif reaction_product == 'p':
                daughter = (df[mass][j]+proton_mass)*MeV
            if model in ('AME20','AME16' ,'I261','I261plusAME20'):
                if pd.isna(df['uncertainty_I261'][i])==True:
                    uncertainty_proj = 'uncertainty_AME20'
                else:
                    uncertainty_proj = 'uncertainty_I261'
                if pd.isna(df['uncertainty_I261'][j])==True:
                    uncertainty_prod = 'uncertainty_AME20'
                else:
                    uncertainty_prod = 'uncertainty_I261'
                target_sigma = df[uncertainty_proj][i]
                remnant_sigma = df[uncertainty_prod][j]
                total_unc=np.sqrt(target_sigma**2+remnant_sigma**2)
            df.loc[[i],[name]]=parent-daughter
            df.loc[[i],[reaction_unc]]=total_unc
        except Exception as er:
            print(er)  
    return (df)


alpha_mass=4.00260325413
neutron_mass=1.00866491588
proton_mass=1.00782503224
MeV=931.49432

#load the experimental data
data=pd.read_excel('/Users/stynikas/Python_codes/Marjut_final_Rh_Ru_Nb_Y_Zr_Mo/Published_masses_I261.xlsx', index_col=None, header=0)
#move to an array representation
Z_I261 = data["Z"].to_numpy()
A_I261 = data["A"].to_numpy()
ME_I261 = data["ME(JYFLTRAP)"].to_numpy()/1000.
sME_I261 = data["dME(JYFLTRAP)"].to_numpy()/1000.
#load the dataframe
df = pd.read_pickle("Database_complete_updated.pkl")
df['mass_excess_I261plusAME20'] = df['mass_excess_AME20']
df['uncertainty_I261plusAME20'] = df['uncertainty_AME20']  
df['mass_I261plusAME20'] = df['mass_AME20'] 
#create the corresponding columns with appropriate headers for the experimental data
reactions= ['a','n','p','g','2n']
models=['I261']
basics = ['mass_excess','uncertainty','mass']
for basic in basics:
    for model in models:
        name='{}_{}'.format(basic,models)
#create channels Q values
df1=[]
models=['I261','I261plusAME20']
for reaction_proj in reactions:
    for reaction_product in reactions:
        for model in models:
            if reaction_proj!=reaction_product:
                name='{}{}_{}'.format(reaction_proj,reaction_product,model)
                reaction_unc='reaction_unc_{}{}_{}'.format(reaction_proj,reaction_product,model)
                df[name]='nan'
                df[reaction_unc]='nan'
for i in range (0,len(Z_I261)):
    df.loc[(df.Z == int(Z_I261[i])) & (df.A == int(A_I261[i])), 'mass_excess_I261'] = ME_I261[i]
    df.loc[(df.Z == int(Z_I261[i])) & (df.A == int(A_I261[i])), 'uncertainty_I261'] = sME_I261[i]
    df.loc[(df.Z == int(Z_I261[i])) & (df.A == int(A_I261[i])), 'mass_excess_I261plusAME20'] = ME_I261[i]
    df.loc[(df.Z == int(Z_I261[i])) & (df.A == int(A_I261[i])), 'uncertainty_I261plusAME20'] = sME_I261[i]    
#print(tabulate(df[['Z','A','mass_excess_AME16','mass_excess_AME20','g2n_AME20','reaction_unc_g2n_AME20','mass_excess_FRDM','mass_excess_FRDM12','mass_excess_HFB_D1m','mass_excess_HFB_Skyrme','mass_excess_HFB_3d']],tablefmt = 'psql'))
#print (tabulate(df.loc[(df.mass_excess_I261 != 'nan'), 'Z']))
df["mass_excess_I261"].astype(float)
df.loc[(df.mass_excess_I261 != 'nan'),'mass_I261']= df["A"].astype(float)+(df["mass_excess_I261"].astype(float)/MeV)
df.loc[(df.mass_excess_I261 != 'nan'),'mass_I261plusAME20']= df["A"].astype(float)+(df["mass_excess_I261"].astype(float)/MeV)
#result_df= df[pd.isna(df['mass_excess_I261'])==False]

for Z_num in range (0,110):
    print(Z_num)
    for A_num in range(0,350):
        #(a,n)
        i=[df.index.values[(df['Z'] == Z_num) &(df['A'] == A_num)]]
        j=[df.index.values[(df['Z'] == (Z_num+2)) & (df['A'] == (A_num+3))]]
        i = np.asarray(i, dtype=np.int32)
        j = np.asarray(j, dtype=np.int32)
        if i.size!=0:
            #(a,n)
            if j.size!=0:
                i=int(i[0])
                j=int(j[0])
                df=calculate_Qvalues('a','n','an',i,j,df)

            #(a,g)
            j=[df.index.values[(df['Z'] == (Z_num+2)) & (df['A'] == (A_num+4))]]
            j = np.asarray(j, dtype=np.int32) 
            if j.size!=0:     
                j=int(j[0])
                df=calculate_Qvalues('a','g','ag',i,j,df)
            
            #(a,n)
            j=[df.index.values[(df['Z'] == (Z_num+2)) & (df['A'] == (A_num+3))]]
            j = np.asarray(j, dtype=np.int32)
            if j.size!=0:  
                j=int(j[0])
                df=calculate_Qvalues('a','n','an',i,j,df)

            #(n,g)
            j=[df.index.values[(df['Z'] == (Z_num)) & (df['A'] == (A_num+1))]]
            j = np.asarray(j, dtype=np.int32)
            if j.size!=0:
                j=int(j[0])
                df=calculate_Qvalues('n','g','ng',i,j,df)

            #(n,p)
            j=[df.index.values[(df['Z'] == (Z_num)+1) & (df['A'] == (A_num))]]
            j = np.asarray(j, dtype=np.int32)
            if j.size!=0:
                j=int(j[0])
                df=calculate_Qvalues('n','p','np',i,j,df)
            #(n,a)
            j=[df.index.values[(df['Z'] == (Z_num)-1) & (df['A'] == (A_num)-3)]]
            j = np.asarray(j, dtype=np.int32)
            if j.size!=0:
                j=int(j[0])
                df=calculate_Qvalues('n','p','np',i,j,df)
            #(p,n)
            j=[df.index.values[(df['Z'] == (Z_num)-1) & (df['A'] == (A_num))]]
            j = np.asarray(j, dtype=np.int32)
            if j.size!=0:
                j=int(j[0])
                df=calculate_Qvalues('p','n','pn',i,j,df)
            #(g,n)=-Sn
            j=[df.index.values[(df['Z'] == (Z_num)) & (df['A'] == (A_num)-1)]]
            j = np.asarray(j, dtype=np.int32)
            if j.size!=0:
                j=int(j[0])
                df=calculate_Qvalues('g','n','gn',i,j,df)
            #(g,2n)=-Sn
            j=[df.index.values[(df['Z'] == (Z_num)) & (df['A'] == (A_num)-2)]]
            j = np.asarray(j, dtype=np.int32)
            if j.size!=0:
                j=int(j[0])
                df=calculate_Qvalues('g','2n','g2n',i,j,df)
result_df= df[pd.isna(df['mass_excess_I261'])==False]
#print (tabulate(result_df[['Z','A','mass_excess_I261','uncertainty_I261','uncertainty_AME20','ng_I261','ng_AME20','reaction_unc_ng_I261','reaction_unc_ng_AME20']]))
result_df = df[df['Z']==45]
print('Z','A','mass_excess_I261','mass_excess_AME20','mass_excess_I261plusAME20','uncertainty_I261','uncertainty_AME20','uncertainty_I261plusAME20','ng_I261','ng_AME20','ng_I261plusAME20','reaction_unc_ng_I261','reaction_unc_ng_AME20','reaction_unc_ng_I261plusAME20')
print (tabulate(result_df[['Z','A','mass_excess_I261','mass_excess_AME20','mass_excess_I261plusAME20','uncertainty_I261','uncertainty_AME20','uncertainty_I261plusAME20','ng_I261','ng_AME20','ng_I261plusAME20','reaction_unc_ng_I261','reaction_unc_ng_AME20','reaction_unc_ng_I261plusAME20']]))

df.to_pickle('Database_complete_plus_I261.pkl')
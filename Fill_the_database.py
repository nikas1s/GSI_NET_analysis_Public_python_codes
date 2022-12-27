import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import fortranformat as ff
from os import listdir
import glob
from tabulate import tabulate
from scipy.optimize import curve_fit

def func(x, a, b, c):
     return a*x**2. + b*x + c



def calculate_Qvalues(reaction_proj,reaction_product,name,i,j,df):
    #constants
    alpha_mass=4.00260325413
    neutron_mass=1.00866491588
    proton_mass=1.00782503224
    MeV=931.49432
    models=['AME20','FRDM','HFB_Skyrme','HFB_3d','HFB_D1m','AME16','FRDM12']
    for model in models:
        name='{}{}_{}'.format(reaction_proj,reaction_product,model)
        mass='mass_{}'.format(model)
        reaction_unc='reaction_unc_{}{}_{}'.format(reaction_proj,reaction_product,model)
        try:
            if reaction_proj == 'a':
                parent = (df[mass][i]+alpha_mass)*MeV
            elif reaction_proj=='n':
                parent = (df[mass][i]+neutron_mass)*MeV
            elif reaction_proj=='p':
                parent = (df[mass][i]+proton_mass)*MeV
            elif reaction_proj=='g':
                parent = (df[mass][i])*MeV

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
            df.loc[[i],[name]]=parent-daughter
        except Exception as er:
            print(er) 
             
        try:   
            if model in ('AME20','AME16'):
                uncertainty = 'uncertainty_{}'.format(model)
                target_sigma = df[uncertainty][i]
                remnant_sigma = df[uncertainty][j]
                total_unc=np.sqrt(target_sigma**2+remnant_sigma**2)
                df.loc[[i],[reaction_unc]]=total_unc
        except Exception as er:
            print(er)

    return (df)


MeV=931.49432
lformat=ff.FortranRecordReader('a1,i3,i5,i5,i5,1x,a3,a4,1x,a13,a11,a11,a9,1x,a2,a11,a9,1x,i3,1x,a12,a11')
file1 = open('/Users/stynikas/Data/mass16.txt', 'r')
count = 0
AME16=[]
while True:
    count += 1
    line = file1.readline()
    if count >39 and count<3437:
        line=lformat.read(line)
        AME16.append(line)
    if not line:
        break
file1.close()
file2 = open('/Users/stynikas/Data/mass20.txt','r')
lformat=ff.FortranRecordReader('a1,i3,i5,i5,i5,1x,a3,a4,1x,a14,a12,a13,1x,a10,1x,a2,a13,a11,1x,i3,1x,a13,a12')
count = 0
AME20=[]
while True:
    count += 1
    line = file2.readline()
    if count >39 and count<3437:
        line=lformat.read(line)
        AME20.append(line)
    if not line:
        break
file2.close()


Z_HFB3D,N_HFB3D,mass_excess_HFB3D=np.loadtxt('HFB_3d.dat',usecols=(0,1,3),unpack=True,skiprows=1)
A_HFB3D=Z_HFB3D+N_HFB3D

Z_FRDM2012,N_FRDM2012,A_FRDM2012 ,mass_excess_FRDM2012=np.loadtxt('FRDM2012.txt',usecols=(0,1,2,14),unpack=True,skiprows=39)


files=glob.glob("/Users/stynikas/Downloads/talys/structure/masses/moller/*")
FRDM=[]
for filename in files:
    file1 = open(filename, 'r')
    for line in file1:
        line=line.strip().split()
        FRDM.append(line)
FRDM = np.array(FRDM)


files=glob.glob("/Users/stynikas/Downloads/talys/structure/masses/hfbd1m/*")
HFBD1m_A=[]
HFBD1m_Z=[]
HFBD1m_ME=[]
for filename in files:
    file1 = open(filename, 'r')
    for line in file1:
        line=line.strip().split()
        HFBD1m_Z.append(int(line[0]))
        HFBD1m_A.append(int(line[1]))
        HFBD1m_ME.append(float(line[3]))
HFBD1m_Z = np.array(HFBD1m_Z)
HFBD1m_A = np.array(HFBD1m_A)
HFBD1m_ME = np.array(HFBD1m_ME)

files=glob.glob("/Users/stynikas/Downloads/talys/structure/masses/hfb/*")
HFB_A=[]
HFB_Z=[]
HFB_ME=[]
for filename in files:
    file1 = open(filename, 'r')
    for line in file1:
        line=line.strip().split()
        line=np.array(line)
        HFB_Z.append(int(line[0]))
        HFB_A.append(int(line[1]))
        HFB_ME.append(float(line[3]))
HFB_Z = np.array(HFB_Z)
HFB_A = np.array(HFB_A)
HFB_ME = np.array(HFB_ME)


#print (HFB[:,1])
print (HFB_Z,HFB_A,HFB_ME[0])

df = pd.DataFrame(data=FRDM,columns=("Z","A","mass","mass_excess_FRDM","dif1","dif2","symbol"))
df=df.drop(columns=["dif1","dif2","mass"])
df=df[["symbol","Z","A","mass_excess_FRDM"]]
df["mass_excess_AME16"]="nan"
df["uncertainty_AME16"]="nan"
df["extrapolated_AME16"]="nan"
df["mass_excess_AME20"]="nan"
df["uncertainty_AME20"]="nan"
df["extrapolated_AME20"]="nan"
df["mass_excess_HFB_3d"]="nan"
df["mass_HFB_3d"]="nan"
df["mass_excess_HFB_Skyrme"]="nan"
df["mass_excess_HFB_D1m"]="nan"
df['mass_excess_FRDM12']="nan"
df=df.sort_values(by=['Z','A'])
s=0
df = df.astype({'Z':'int'})
df = df.astype({'A':'int'})
for i in range (0,len(Z_HFB3D)):
    df.loc[(df.Z == int(Z_HFB3D[i])) & (df.A == int(A_HFB3D[i])), 'mass_excess_HFB_3d'] = mass_excess_HFB3D[i]

for i in range(0,len(mass_excess_FRDM2012)):
    df.loc[(df.Z == int(Z_FRDM2012[i])) & (df.A == int(A_FRDM2012[i])), 'mass_excess_FRDM12'] = mass_excess_FRDM2012[i]

for i in range(0,len(HFB_ME)):
    df.loc[(df.Z == int(HFB_Z[i])) & (df.A == int(HFB_A[i])), 'mass_excess_HFB_Skyrme'] = HFB_ME[i]

for i in range(0,len(HFBD1m_ME)):
    df.loc[(df.Z == int(HFBD1m_Z[i])) & (df.A == int(HFBD1m_A[i])), 'mass_excess_HFB_D1m'] = HFBD1m_ME[i]
    


for i in range (0,len(AME16)):
    s=s+1
    position=df.index[(df.Z == int(AME16[i][3])) & (df.A == int(AME16[i][4]))]
    if "#" in AME16[i][7]:
        AME16[i][7]=AME16[i][7].replace('#', '')
        AME16[i][8]=AME16[i][8].replace('#', '')
        df.loc[(df.Z == int(AME16[i][3])) & (df.A == int(AME16[i][4])), "extrapolated_AME16"]="yes"
    else:
        df.loc[(df.Z == int(AME16[i][3])) & (df.A == int(AME16[i][4])), "extrapolated_AME16"]="no"
    AME16[i][8]=float(AME16[i][8])/1000.
    AME16[i][7]=float(AME16[i][7])/1000.
    #print(AME16[i][3],AME16[i][4],AME16[i][7],position)
    df.loc[(df.Z == int(AME16[i][3])) & (df.A == int(AME16[i][4])), "mass_excess_AME16" ]=AME16[i][7]
    df.loc[(df.Z == int(AME16[i][3])) & (df.A == int(AME16[i][4])), "uncertainty_AME16"]=AME16[i][8]

for i in range (0,len(AME20)):
    s=s+1
    position=df.index[(df.Z == int(AME20[i][3])) & (df.A == int(AME20[i][4]))]
    if "#" in AME20[i][7]:
        AME20[i][7]=AME20[i][7].replace('#', '')
        AME20[i][8]=AME20[i][8].replace('#', '')
        df.loc[(df.Z == int(AME20[i][3])) & (df.A == int(AME20[i][4])), "extrapolated_AME20"]="yes"
    else:
        df.loc[(df.Z == int(AME20[i][3])) & (df.A == int(AME20[i][4])), "extrapolated_AME20"]="no"
    AME20[i][8]=float(AME20[i][8])/1000.
    AME20[i][7]=float(AME20[i][7])/1000.
    df.loc[(df.Z == int(AME20[i][3])) & (df.A == int(AME20[i][4])), "mass_excess_AME20" ]=AME20[i][7]
    df.loc[(df.Z == int(AME20[i][3])) & (df.A == int(AME20[i][4])), "uncertainty_AME20"]=AME20[i][8]


df["neutron_rich"]="nan"
df.loc[(df.mass_excess_HFB_D1m != 'nan'),'mass_HFB_D1m']= df["A"].astype(float)+(df["mass_excess_HFB_D1m"].astype(float)/MeV)
df.loc[(df.mass_excess_HFB_Skyrme != 'nan'),'mass_HFB_Skyrme']= df["A"].astype(float)+(df["mass_excess_HFB_Skyrme"].astype(float)/MeV)
df.loc[(df.mass_excess_HFB_3d != 'nan'),'mass_HFB_3d']= df["A"].astype(float)+(df["mass_excess_HFB_3d"].astype(float)/MeV)
df.loc[(df.mass_excess_AME20 != 'nan'),'mass_AME20']= df["A"].astype(float)+(df["mass_excess_AME20"].astype(float)/MeV)
df.loc[(df.mass_excess_AME16 != 'nan'),'mass_AME16']= df["A"].astype(float)+(df["mass_excess_AME16"].astype(float)/MeV)
df.loc[(df.mass_excess_FRDM != 'nan'),'mass_FRDM']= df["A"].astype(float)+(df["mass_excess_FRDM"].astype(float)/MeV)
df.loc[(df.mass_excess_FRDM12 != 'nan'),'mass_FRDM12']= df["A"].astype(float)+(df["mass_excess_FRDM12"].astype(float)/MeV)
stable=np.loadtxt("/Users/stynikas/Data/Stable_Nuclides.txt",unpack=True)
popt, pcov = curve_fit(func, stable[0], stable[1])
df["position_help"]=func(df["A"].astype(int),popt[0],popt[1],popt[2])
df.loc[df["position_help"]<df["Z"].astype(float),"neutron_rich"]="p-rich"
df.loc[df["position_help"]>df["Z"].astype(float),"neutron_rich"]="n-rich"
df["Z"]=df["Z"].astype(int)
df["A"]=df["A"].astype(int)
for i in range (0,len(stable[1])):
    position=df.index[df.Z.isin([stable[1][i]]) & df.A.isin([stable[0][i]])]
    df["neutron_rich"][position]="stable"

input_data=df[((df['neutron_rich'] == "n-rich") & (df['extrapolated_AME16'] == "no")) | ((df['neutron_rich'] == "n-rich") & (df['extrapolated_AME16'] == "nan")) | ((df['neutron_rich'] == "n-rich") & (df['extrapolated_AME16'] == "yes"))]
        
        

reactions= ['a','n','p','g','2n']
models=['AME20','AME16','FRDM','FRDM12','HFB_Skyrme','HFB_3d','HFB_D1m']
for reaction_proj in reactions:
    for reaction_product in reactions:
        for model in models:
            if reaction_proj!=reaction_product:
                name='{}{}_{}'.format(reaction_proj,reaction_product,model)
                df[name]='nan'
            if model=='AME20' or model=='AME16':    
                reaction_unc='reaction_unc_{}{}_{}'.format(reaction_proj,reaction_product,model)
                df[reaction_unc]='nan'
df["mass_FRDM"].astype(float)
df["mass_FRDM12"].astype(float)
df["mass_AME20"].astype(float)

df["mass_HFB_3d"].astype(float)
df["Z"]=df["Z"].astype(int)
df["A"]=df["A"].astype(int)
for Z_num in range (0,110):
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

    print(Z_num)

df.to_pickle('Database_complete_updated.pkl')    #to save the dataframe, df to 123.pkl
#print(tabulate(df[['Z', 'A','mass_excess_AME20','mass_excess_FRDM12','ag_FRDM12','g2n_HFB_3d','g2n_FRDM','g2n_FRDM12','mass_excess_HFB_3d','mass_HFB_3d','an_HFB_3d','ag_HFB_3d','mass_excess_FRDM','mass_excess_AME20']],tablefmt = 'psql'))
print(tabulate(df[['Z','A','mass_excess_AME16','mass_excess_AME20','g2n_AME20','reaction_unc_g2n_AME20','mass_excess_FRDM','mass_excess_FRDM12','mass_excess_HFB_D1m','mass_excess_HFB_Skyrme','mass_excess_HFB_3d']],tablefmt = 'psql'))



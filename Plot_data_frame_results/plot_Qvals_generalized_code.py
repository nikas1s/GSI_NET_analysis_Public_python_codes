import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tabulate import tabulate

#choose the mass number range of nuclei to be plotted
A_lims=[[107,123],[112,120],[100,110],[100,108],[103,109],[108,115]]

#color blind friendly pallete from http://colorbrewer2.org/
colors=['#1f78b4','#b2df8a','#33a02c','#fb9a99','#fdbf6f','#e31a1c','#a6cee3']


#Load the pandas dataframe corresponding to data from I261 experiment
df = pd.read_pickle('/Users/stynikas/Python_codes/databases_codes/Database_complete_plus_I261.pkl')
#position counter set initial to -1 update in beginning of the loop counter+=1
counter=-1


plt.rc('figure', figsize=(5., 7))
for Z in [45,44,41,39,40,42]:
    #define the models to be plotted
    models=['FRDM12','I261','AME20','HFB_Skyrme','HFB_3d','I261plusAME20']
    #projectile and product of the reaction 'g' 'n' will i.e. plot the g,n Q - values
    # you can choose between n,2n,g,p,a keep in mind (a,a) (elements in the diagonal)
    # do not exist since Qvalue is 0 :)
    reaction_proj='g'
    reaction_product='2n'
    #create a figure with 3 subplots
    ax1 = plt.subplot(3, 1,2)
    ax0 = plt.subplot(3, 1,1, sharex=ax1)
    ax2 = plt.subplot(3, 1,3, sharex=ax1)
    counter+=1
    mod_num=-1
    for model in models:
        #generate the names from the database
        mod_num+=1
        name='{}{}_{}'.format(reaction_proj,reaction_product,model)
        reaction_unc='reaction_unc_{}{}_{}'.format(reaction_proj,reaction_product,model)
        mass='mass_{}'.format(model)
        # Create an instance of the dataframe to manipulate called r_df
        r_df=df
        #choose only the elements where atomic number Z = df[Z]
        result_df= r_df[df['Z']==Z]
        result_df=result_df.sort_values(by=['A'])
        A=result_df['A'].to_numpy()
        print (A)
        Qv=result_df[name].to_numpy()
        print (Qv)
        if mod_num==1:
            r_df=df[pd.isna(df[mass])==False]
            result_df= r_df[df['Z']==Z]
            result_df=result_df.sort_values(by=['A'])
            A_exp=result_df['A'].to_numpy()
            Qv_exp=result_df[name].to_numpy()

        if model in ('I261','AME20','AME16','I261plusAME20'):
            sigma_Qv=result_df[reaction_unc].to_numpy()
            sigma_Qv[sigma_Qv=='nan']=0
        if model in ('FRDM12'):
            stad=Qv
            stad_A=A
            stad[stad=='nan']=0
        if model in ('AME20'):
            stad_s=sigma_Qv
        Qv[Qv=='nan']=0
        #plot error bars for all the models that have errror bars (models 2 & 5)
        if mod_num==2 or mod_num==5:
            ax0.errorbar(A,Qv,yerr=sigma_Qv,xerr=0)
            if mod_num==2:
                ax1.errorbar(A,Qv-stad,yerr=sigma_Qv,xerr=0,capsize=5,capthick=2,color=colors[mod_num])
                ax0.errorbar(A,Qv,yerr=sigma_Qv,xerr=0, capsize=5,capthick=2,label='{}'.format(model),color=colors[mod_num])
            else:
                ax1.errorbar(A, Qv-stad, yerr = sigma_Qv, xerr = 0,capsize=2,capthick=1,color=colors[mod_num])
                ax0.errorbar(A,Qv,yerr=sigma_Qv,xerr=0,capsize=2,capthick=1,label='{}'.format(model),color=colors[mod_num])
                ax2.scatter(A,(stad_s-sigma_Qv)/(stad_s+0.000000001)*100,color=colors[mod_num])
        #else plot single lines with no error bars
        elif mod_num==0 or mod_num==3 or mod_num==4:
            ax0.plot(A,Qv,label='{}'.format(model),color=colors[mod_num])
            ax1.plot(A,Qv-stad,color=colors[mod_num])

        ax0.set_xlim(A_lims[counter][0],A_lims[counter][1])
        ax0.set_ylabel('MeV')
        ax1.set_ylabel('$ME-ME_{FRDM}$')
        ax2.set_xlabel('A (mass number)')
        ax2.set_ylabel('relative $\sigma$ reduction')
        ax1.set_ylim(-1.5,1.5)
    ax0.legend(loc='lower left')
    #ax0.set_ylim(-5,10)
    ax0.scatter(A_exp, Qv_exp, s=100, facecolors='none', edgecolors='r')
    #plt.show()
    #
    ax0.set_title('Qvalue for ({}{}) for Z={}'.format(reaction_proj,reaction_product,Z))
    plt.savefig('Plots_examples/{}{}_{}.pdf'.format(reaction_proj,reaction_product,Z))
    plt.close()
    #plt.show()

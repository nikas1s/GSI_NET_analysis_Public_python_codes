# This routine Fits in a 7 parameter formula (a0,a6) n,g reaction rates (how fast a reaction happens)
# Following the reaclib conventions more info  https://reaclib.jinaweb.org/docs/reaclibFormat.pdf
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pylab as py
import random
import csv
from itertools import islice
import reaclib_format as rf
import os

# list of elements needed for the symbols#
elems=(['nt','h','he','li', 'be', ' b', ' c', ' n', ' o', ' f','ne', 'na', 'mg', 'al', 'si', ' p', ' s','cl', 'ar', ' k', 'ca', 'sc', 'ti', ' v','cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn','ga', 'ge', 'as', 'se', 'br', 'kr', 'rb','sr', ' y', 'zr', 'nb', 'mo', 'tc', 'ru','rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb','te', ' i', 'xe', 'cs', 'ba', 'la', 'ce','pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb','dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf','ta', ' w', 're', 'os', 'ir', 'pt', 'au','hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn','fr', 'ra', 'ac', 'th', 'pa', ' u', 'np','pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm','md', 'no', 'lr', 'rf', 'db', 'sg', 'bh','hs', 'mt','ds','rg','cn','nh','fl','mc','116','ts','og' ])

#============== define the function based on the reaclib format==========================
def fitFunc(T,a0,a1,a2,a3,a4,a5,a6):

    return (a0+a1*T**(-1)+a2*T**(-1./3.)+a3*T**(1./3.)+a4*T+a5*T**(5./3.)+a6*np.log10(T))
#==============================================================================

#==============================================================================
#==================Main Programm===============================================
#==============================================================================


#===========Define a list of the folder names to be accesed====================
path=['rates_example/']
labels=['AME20']
Q_values=[] #list to input Q values
spin=[] #list to input the spins
#optional give the Z and the A range for the nuclei
A1=[[80,83]]
Z1=[30]
for counter in (0,len(path)):
    path1="fittings/{}".format(labels[counter])
    os.mkdir(path1)
    for i in range (0,len(Z1)):
        for j in range (A1[i][0],A1[i][1]):
            Z=Z1[i] # atomic number
            A=j # mass number
            N=A-Z #neutron number
            Aplusone=A+1 #neghboring nucleus (A+1)
            a=[]
            ll=0
            ############read the file containing the reactions########################
            filename=("{}{}_{}#{}.bin".format(path[counter],A,Z,labels[counter]))
            try:
                #try to load the T and reaction rate as xdata and ydata
                xdata,ydata=py.loadtxt(filename,skiprows=2,usecols=(0,2),unpack=True)
            except:
                print ('did not find the file {}'.format(filename))
            try:
                #try to load the file that contains spin information
                filename2=("{}{}_{}#{}.out".format(path[counter],A,Z,labels[counter]))
                file=open(filename2,'r')
            except:
                print ('did not find the out file {}'.format(filename2))
            ###########################################################################

            #patterns in the file to identify and read the corresponding needed data
            pattern = 'Q(n,g'
            pattern2='   0     0.0000  '
            pattern3=' massnucleus {}'.format(Z)
            ##########################################################################

            del Q_values[:]
            del spin[:]
            #loop over all the lines of the file
            for line in file:
                #Find the Q value
                if pattern in line:
                    print ("I MADE IT read q value")
                    Q_values.append(line.strip().split(':')[1])
                    Q_values[0]=float(Q_values[0])
                #find the spin
                if pattern2 in line:
                    spin.append(line.split(' ')[10])
                    spin[0]=float(spin[0])
            try:
                #Read the neighboring nuclei data
                filename_A_plus_one=("{}{}_{}#{}.out".format(path[counter],Aplusone,Z,labels[counter]))
                file2 = open(filename_A_plus_one, "r")
            except:
                print ('did not find a plus one {}'.format(filename_A_plus_one))
            
            for line in file2:
                #find the spin
                 if pattern2 in line:

                    spin.append(line.split(' ')[10])
                    spin[2]=float(spin[2])
            try:
                #move to log representation to account for the multiple orders of magnitude
                ydata=np.log(ydata)
                print ('here')
                #initialize the factors a as a random resonable number
                del a[:]
                a.append(np.random.uniform(-1500.,1500.))
                a.append(np.random.uniform(-1500.,1500.))
                a.append(np.random.uniform(-1500.,1500.))
                a.append(np.random.uniform(-1500.,1500.))
                a.append(np.random.uniform(-1500.,1500.))
                a.append(np.random.uniform(-1500.,1500.))
                a.append(np.random.uniform(-1500.,1500.))
                fitParams1=a
                for l in range (1,10000):

                        #=========== Try to fit the data if error occurs the process continues with a new set of data==============
                        fitParams, fitCovariance = curve_fit(fitFunc, xdata,ydata,[fitParams1[0],fitParams1[1],fitParams1[2],fitParams1[3],fitParams1[4],fitParams1[5],fitParams1[6]])
                        #=======Initialize the sum and choose the minimum chi2 from the ones provided from the initial guesses====
                        sum_ch=0.0
                        for k in range (0,30):
                         chi2=(((fitFunc(xdata[k], *fitParams))-(ydata[k]))**2)/(ydata[k])
                         sum_ch =sum_ch + chi2
                        if l==1:
                            sum_ch_compare=1.7976931348623157e+308
                        if sum_ch<sum_ch_compare:
                            sum_ch_compare=sum_ch
                            fitParams1=fitParams
                            print('I am reducing')
            except:
                print("some value is 0")
            print (sum_ch)
            #==========print the parameters plot the results========================
            #print fitParams1
            plt.clf()
            r1=rf.reaclib_entry(4,[fitParams1[0],fitParams1[1],fitParams[2],fitParams1[3],fitParams1[4],fitParams1[5],fitParams1[6]],["n","{}{}".format(elems[Z],A),"{}{}".format(elems[Z],Aplusone),"=====","=====","====="],"SNPD",False,False,False,Q_values[0])
            if counter!=2:
                r1.print_to_file("fittings/{}/fittings_{}_{}".format(labels[counter],Z,N))
            if counter==2:
                r1.print_to_file("fittings/{}/fittings_{}_{}".format(labels[4],Z,N))
            #=======================================================================
            #=============Calculate reverse rate====================================
            #=======================================================================
            A_f=float(A)
            coeffs_rev=[]
            del coeffs_rev[:]
            c=(A_f/(A_f+1.))**(1.5)
            a=((2.*(spin[0]+1.))/(spin[2]+1.))
            b=3./(2.*(spin[2]+1.))*9.8685e9
            balance_factor=a*b*c

            for ii in range (0,7):
                 coeffs_rev.append(fitParams1[ii])
            coeffs_rev[0]=coeffs_rev[0]+np.log10(balance_factor)
            coeffs_rev[1]=coeffs_rev[1]+11.6045*Q_values[0]*(-1.)
            coeffs_rev[6]=coeffs_rev[6]+1.5
            r1=rf.reaclib_entry(2,[coeffs_rev[0],coeffs_rev[1],coeffs_rev[2],coeffs_rev[3],coeffs_rev[4],coeffs_rev[5],coeffs_rev[6]],["{}{}".format(elems[Z],Aplusone),"n","{}{}".format(elems[Z],A),"=====","=====","====="],"SNPD",False,False,True,-Q_values[0])

            r1.print_to_file("fittings/{}/fittings_{}_{}".format(labels[counter],Z,N))

            #Plot the rest of the stuff
            print ("so far so good")
            plt.xlabel('Temperature (GK)', fontsize = 16)
            #===========================================================================
            plt.ylabel('RR', fontsize = 16)
            xdata1=np.linspace(0.001,10,1000)
            plt.semilogy(xdata1,py.exp(fitFunc(xdata1, *fitParams1)),linewidth=10.,alpha=0.15)
            plt.semilogy(xdata1[20:],py.exp(fitFunc(xdata1[20:], *coeffs_rev)),linewidth=6.,alpha=0.20)
            plt.semilogy(xdata, py.exp(ydata), 'ro')
            #plt.ylim(min(py.exp(fitFunc(xdata1, *fitParams1))),max(py.exp(fitFunc(xdata1, *fitParams1))))
            plt.savefig("plots/rlib16AGMaximes_fit_{}_{}_{}.png".format(labels[counter],Z,N))
      #except:
        #print("SOME ERROR")

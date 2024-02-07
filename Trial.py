import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import os
import math
from veldists import VelDist
from detectors.nai import DAMA
from constants import *

new_directory = r"C:\Users\SAMARTH\OneDrive\Desktop\Python Code\sgm"

# Change the current working directory
os.chdir(new_directory)

Det = DAMA()



mX = 690 * keV


def Error(ER):
    k = math.erf((ER+26)*(1/13.7))
    return k 

def DMRate(mX,ER,T,sig,veldist):
    Beta = 10**(-3)
    Phi_DM = 7.8*10**(-4)*(1/mX)*Beta
    g = 1
    sig_PE = mX**2
    sigE = np.sqrt((0.0256+0.0003256*ER))
    sig_AE = sig_PE*(g**2)/Beta *((3*mX**2)/16*np.pi*511*511*0.007297)*(1-Beta**(2/3)/3)
    MT = 365 #Exposure in days

    dsigdER = Phi_DM*sig_AE*Error(ER)*(1/np.sqrt(2*np.pi))*sigE*np.exp(-(ER-mX)**2/2*(sigE**2))*MT

    return dsigdER

ER = np.arange(10,100,0.1) # observed energy, units of keV_ee
plt.plot(ER,[Det.dRdE(DMRate,e)for e in ER])
plt.ylabel(r"Observation rate [cpd/kg/keV$_{\rm ee}$]")
plt.xlabel(r"Energy [keV$_{\rm ee}$]")
#plt.xlim(0,16)
#plt.ylim(0,0.026)
plt.grid()
plt.show()
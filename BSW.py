from dmmodel import DMModel
import numpy as np 
import math
import os 
from constants import *

new_directory = r"C:\Users\SAMARTH\OneDrive\Desktop\Python Code\sgm"

# Change the current working directory
os.chdir(new_directory)

class SIWIMPBSW(DMModel):
    def Error(self,ER):
        return math.erf(ER+26)*(1/13.7)

    def dRdER(self,mX, ER):
        Beta = 10 ** (-3)
        Phi_DM = 7.8*10**(-4)*(1/mX)*Beta
        sig_PE = mX**2
        g = 1
        sigE = np.sqrt((0.0256+0.0003256*ER))
        sig_AE = sig_PE*(g**2)/Beta *((3*mX**2)/16*np.pi*511*511*0.007297)*(1-Beta**(2/3)/3)
        MT = 365 #Exposure in days

        dsigdER = Phi_DM*sig_AE*SIWIMPBSW.Error(ER)*(1/np.sqrt(2*np.pi))*sigE*np.exp(-(ER-mX)**2/2*(sigE**2))*MT

        return dsigdER
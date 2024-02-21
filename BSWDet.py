import numpy as np
from detector import Detector
from targets.na import Na
from targets.i import I
from constants import *
from BSW import SIWIMPBSW



class DAMA_BSW(Detector):
    def Nuclei(self):
        return[Na(),I()]
    
    def ER_E(self,E):
        #Quenching Factor

        return [E*keV/0.3,E*keV/0.09]
    
    def dERdE(self, E):
        return [1/0.3,1/0.09]
    
    def ROI(self):
        return [10,1000]
    
    def Emax(self):
        return 700
    
    def DeltaE(self,E):
        return (0.488*pow(E,0.5))+(0.0091*E)
    
    def Res(self,E1,E2):
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        R = A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
        return R
    def dRdE(self,E,sig,mX):
        
        return lambda E2:SIWIMPBSW.dRdER_a(mX, sig)*self.Res(E,E2)

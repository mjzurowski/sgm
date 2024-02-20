import numpy as np
from detector import Detector
from targets.na import Na
from targets.i import I
from constants import *

class DAMA(Detector):
    def Nuclei(self):
        return [Na(),I()]
    
    def ER_E(self,E):
        ### This function will ultimately need some way to distinguish between NR and ER
        ### For now can just avoid that by creating a separate background spec
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        return [E*keV/0.3, E*keV/0.09]
    
    def dERdE(self,E):
        return [1/0.3, 1/0.09] ### ideally this should just be computed automatically in detector.py
    
    def ROI(self):
        return [10,1000]
    
    def Emax(self):
        return 1000
    
    def DeltaE(self,E):
        return (0.488*pow(E,0.5))+(0.0091*E)
    
    def Res(self,E1,E2):
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        return np.where(E>8,1,0.0429*E+0.657 ) 
    
    def dRdER(self, Target, ER, mX, sig, VelDist):
        if ER!=mX:
            return 0
        else:
            C_a = 1.2 * 10**(19)
            g_ae = sig
            dsigdER = (C_a/(Target.A))*(g_ae**2)*mX*(Target.sigma_PE(ER))
            return dsigdER
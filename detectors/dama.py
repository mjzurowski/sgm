import numpy as np
from detector import Detector
from targets.na import Na
from targets.i import I
from constants import *

class DAMA(Detector):
    def Nuclei(self):
        return [Na(),I()]
    
    def ER_E(self,E):
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
        return 690
    
    def DeltaE(self,E):
        return (0.488*pow(E,0.5))+(0.0091*E)
    
    def Res(self,E1,E2):
        # We assume E1 is the observed energy (E' in accompanying documentation) and E2 is the energy that will be integrated over (E_ee in accompanying documentation)
        A = 1./(np.sqrt(2.*np.pi))
        return A*np.exp(-0.5*pow((E1 - E2), 2.))
    
    def Eff(self,E):
        return np.where(E>8,1,0.0429*E+0.657 ) 
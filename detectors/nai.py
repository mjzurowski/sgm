import numpy as np
from detector import Detector
from targets.na import Na

class DAMA(Detector):
    def Nuclei(self):
        return [Na()]
    
    def ER_E(self,E):
        ### This function will ultimately need some way to distinguish between NR and ER
        ### For now can just avoid that by creating a separate background spec
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        return [E*1E3/0.3]
    
    def dERdE(self,E):
        return [1/0.3]
    
    def ROI(self):
        return [0,6]
    
    def Emax(self):
        return 20
    
    def DeltaE(self,E):
        return 0.4427*pow(E,0.5)
    
    def Res(self,E1,E2):
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        return np.where(E>8,1,0.0429*E+0.657 ) 
import numpy as np
from detector import Detector
from targets.na import Na
from targets.i import I
from constants import *
from detectors.dama import DAMA


class DAMA_BSW(Detector):
    
    def DeltaE(self,E):
        sigma_0 = 0.16
        return np.sqrt(sigma_0**(2)+2.96*0.11*E)
    
    def Res(self,E1,E2):
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        R = A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
        return R
    

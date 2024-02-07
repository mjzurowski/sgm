import numpy as np
from detector import Detector
from targets.na import Na
from targets.i import I
from constants import *

joo_qf = np.loadtxt("./detectors/na_qfs/joo.dat")  ## list of E_NR as a function of E_ee from Joo QF measurements

class COSINE(Detector):
    def Nuclei(self):
        return [Na(),I()]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV based on Joo QFs
        """
        return [np.interp(E,joo_qf[:,0],joo_qf[:,1])*keV, E*keV/0.05]
    
    def dERdE(self,E):
        return [np.interp(E,joo_qf[:,0],joo_qf[:,1])/E, 1/0.05] ### ideally this should just be computed automatically in detector.py
    
    def ROI(self):
        return [0,6]
    
    def Emax(self):
        return 20
    
    def DeltaE(self,E):
        return (0.3171*pow(E,0.5))+(0.008189*E)
    
    def Res(self,E1,E2):
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        return np.where(E>2,1,0.2*E+0.6) 
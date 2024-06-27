import numpy as np
from detector import Detector
from targets.ge import Ge
from constants import *

class SCDMS_GeHV(Detector):
    def __init__(self, volt, shell_model="Fitz"):
        ## allow for initialisation with different shell models
        self.shell_model = shell_model
        self.V = volt # voltage detector is run at in volts

    def Nuclei(self):
        return [Ge(self.shell_model)]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        eps = 3.0 #eV
        Y = 0.2 # ionisation yield - need to think about how to get a nicer expression for this
        scale = 1+Y*e*self.V/eps
        return [E*keV/scale] # better yet, this could/should be driven by calibration method, or some kind of interpolation
    
    def dERdE(self,E):
        eps = 3.0 #eV
        Y = 0.2 # ionisation yield - need to think about how to get a nicer expression for this
        scale = 1+Y*e*self.V/eps
        return [1/scale] ### ideally this should just be computed automatically in detector.py
    
    def ROI(self):
        return [0,10]
    
    def Emax(self):
        return 20
    
    def DeltaE(self,E):
        A = 5E-3
        B = 0.7 
        sig_E = 10
        return np.sqrt(B*E/keV+pow(A*E/keV,2)+pow(sig_E/keV,2))*keV # these resolutions are defined for eV so need to convert to eV for comp of DeltaE
    
    def Res(self,E1,E2):
        # We assume E1 is the observed energy (E' in accompanying documentation) and E2 is the energy that will be integrated over (E_ee in accompanying documentation)
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        return np.where(E>0.15,0.8,0.8*E/0.15) 
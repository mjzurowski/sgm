import numpy as np
from detector import Detector
from targets.si import *
from constants import *

"""
Define a "detector" object to combine all the isotopes of Ge without needing a particular detector
"""

class Si(Detector):
    def __init__(self, volt, shell_model="Fitz"):
        ## allow for initialisation with different shell models
        self.si28 = Si28(shell_model)
        self.si29 = Si29(shell_model)
        self.si30 = Si30(shell_model)

    def Nuclei(self):
        return [[self.si28,0.922],[self.si29,0.047],[self.si30,0.031]]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        return [E*keV,E*keV,E*keV] # this could alternatively be driven by calibration process
    
    def dERdE(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output derivative of ER wrt E_obs, units [keV]/[keV_0]
        """
        return np.ones(len(self.Nuclei())) # interpolate derivative at the value we want 
    
    def ROI(self):
        return [0,100]
    
    def Emax(self):
        return 200
    
    def DeltaE(self,E):
        # Not used
        return 0 # these resolutions are defined for eV so need to convert to keV for comp of DeltaE
    
    def Res(self,E1,E2):
        # Assume perfect resolution
        np.where(E1==E2,1,0)
    
    def Eff(self,E):
        # Assume perfect efficiency
        return 1

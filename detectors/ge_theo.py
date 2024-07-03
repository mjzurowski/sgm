import numpy as np
from detector import Detector
from targets.ge import Ge70, Ge72, Ge73, Ge74, Ge76
from constants import *

"""
Define a "detector" object to combine all the isotopes of Ge without needing a particular detector
"""

class Ge(Detector):
    def __init__(self, volt, shell_model="Fitz"):
        ## allow for initialisation with different shell models
        self.ge70 = Ge70(shell_model)
        self.ge72 = Ge72(shell_model)
        self.ge73 = Ge73(shell_model)
        self.ge74 = Ge74(shell_model)
        self.ge76 = Ge76(shell_model)

    def Nuclei(self):
        return [self.ge70, self.ge72, self.ge73, self.ge74, self.ge76]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        return [E*keV,E*keV,E*keV,E*keV,E*keV] # this could alternatively be driven by calibration process
    
    def dERdE(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output derivative of ER wrt E_obs, units [keV]/[keV_0]
        """
        # Here we account for their abundance fraction. We want to scale kg of isotope per kg of material
        A_av = 0.205*self.ge70.A()+0.274*self.ge72.A()+0.0776*self.ge73.A()+0.365*self.ge74.A()+0.0775*self.ge76.A() # get the average mass number
        return [0.205*self.ge70.A()/A_av,0.274*self.ge72.A()/A_av,0.0776*self.ge73.A()/A_av,0.365*self.ge74.A()/A_av,0.0775*self.ge76.A()/A_av] # interpolate derivative at the value we want 
    
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

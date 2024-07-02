import numpy as np
from detector import Detector
from targets.ge import Ge
from targets.si import Si
from constants import *

def Y_Si(ER):
    """
    Ionisation yield for Si. Taken from 10.1103/PhysRevLett.131.091801
    [ER] = [eV]
    """
    Y_10keV = 0.302
    B = 0.261
    return Y_10keV*pow(ER*1E-4,B)

def Y_Ge(ER):
    """
    Ionisation yield for Ge. Taken from the Lindhard assumptions given in arxiv 1304.6773
    """
    Z = 32
    A = 72.64
    k = 0.133*pow(Z,2/3)*pow(A,-0.5)
    ep = 11.5*ER*pow(Z,-7/3)
    g = 3*pow(ep,0.15)+0.7*pow(ep,0.6)+ep
    return k*g/(1+k*g)

class GeHV(Detector):
    def __init__(self, volt, shell_model="Fitz"):
        ## allow for initialisation with different shell models
        self.shell_model = shell_model
        self.V = volt # voltage detector is run at in volts

        # construct arrays for E_obs --> ER interpolation
        eps = 3.0 #eV
        self.ER_samp = np.arange(0,100*keV,100) # recoil energy in eV range up to 100 keV in steps of 100 eV
        self.E_samp = (Y_Ge(self.ER_samp)*volt/eps)*self.ER_samp + self.ER_samp # observed energy in eV

    def Nuclei(self):
        return [Ge(self.shell_model)]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        ER = np.interp(E*keV,self.E_samp,self.ER_samp)
        return [ER] # this could alternatively be driven by calibration process
    
    def dERdE(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output derivative of ER wrt E_obs, units [keV]/[keV_0]
        """
        deriv = np.gradient(self.ER_samp,self.E_samp) # get the derivative of ER vs E_obs
        return [np.interp(E*keV,self.E_samp,deriv)] # interpolate derivative at the value we want 
    
    def ROI(self):
        return [0,10]
    
    def Emax(self):
        return 20
    
    def DeltaE(self,E):
        # Lets assume CDMSlite values: arxiv 1911.11905. These are for a Ge iZIP, so we'll at least adjust the baseline res to what we hope for
        A = 5E-3
        B = 0.7 
        sig_E = 10
        return np.sqrt(B*E/keV+pow(A*E/keV,2)+pow(sig_E/keV,2)) # these resolutions are defined for eV so need to convert to keV for comp of DeltaE
    
    def Res(self,E1,E2):
        # We assume E1 is the observed energy (E' in accompanying documentation) and E2 is the energy that will be integrated over (E_ee in accompanying documentation)
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        return np.where(E>0.15,0.8,0.8*E/0.15) 
    
class SiHV(Detector):
    def __init__(self, volt, shell_model="Fitz"):
        ## allow for initialisation with different shell models
        self.shell_model = shell_model
        self.V = volt # voltage detector is run at in volts

        # construct arrays for E_osb --> ER interpolation
        eps = 3.82 #eV
        self.ER_samp = np.arange(0,100*keV,100) # recoil energy in eV range up to 100 keV in steps of 100 eV
        self.E_samp = (Y_Si(self.ER_samp)*volt/eps)*self.ER_samp + self.ER_samp # observed energy in eV

    def Nuclei(self):
        return [Si(self.shell_model)]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        ER = np.interp(E*keV,self.E_samp,self.ER_samp)
        return [ER] # this could alternatively be driven by calibration process
    
    def dERdE(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output derivative of ER wrt E_obs, units [keV]/[keV_0]
        """
        deriv = np.gradient(self.ER_samp,self.E_samp) # get the derivative of ER vs E_obs
        return [np.interp(E*keV,self.E_samp,deriv)] # interpolate derivative at the value we want 
    
    def ROI(self):
        return [0,10]
    
    def Emax(self):
        return 20
    
    def DeltaE(self,E):
        # Lets assume CDMSlite values: arxiv 1911.11905. These are for a Ge iZIP, so we'll at least adjust the baseline res to what we hope for
        A = 5E-3
        B = 0.7 
        sig_E = 5
        return np.sqrt(B*E/keV+pow(A*E/keV,2)+pow(sig_E/keV,2)) # these resolutions are defined for eV so need to convert to keV for comp of DeltaE
    
    def Res(self,E1,E2):
        # We assume E1 is the observed energy (E' in accompanying documentation) and E2 is the energy that will be integrated over (E_ee in accompanying documentation)
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        # Taken from eyeballing the CDMS Soudan data - should be improved
        return np.where(E>0.15,0.8,0.8*E/0.15) 
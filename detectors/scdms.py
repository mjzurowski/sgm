import numpy as np
from detector import Detector
from targets.ge import *
from targets.si import *
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

        # lets now define the germanium nucleus based on the isotopes
        self.ge70 = Ge70(shell_model)
        self.ge72 = Ge72(shell_model)
        self.ge73 = Ge73(shell_model)
        self.ge74 = Ge74(shell_model)
        self.ge76 = Ge76(shell_model)

        self.V = volt # voltage detector is run at in volts

        # construct arrays for E_obs --> ER interpolation
        eps = 3.0 #eV
        self.ER_samp = np.arange(0,100*keV,100) # recoil energy in eV range up to 100 keV in steps of 100 eV
        self.E_samp = (Y_Ge(self.ER_samp)*volt/eps)*self.ER_samp + self.ER_samp # observed energy in eV

    def Nuclei(self):
        """
        Each entry should be a target nucleus and its abundance in moles per mole of full target
        """
        return [[self.ge70,0.205], [self.ge72,0.274], [self.ge73,0.0776], [self.ge74,0.365], [self.ge76,0.0775]]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        ER = np.interp(E*keV,self.E_samp,self.ER_samp) # interpolate the recoil energy based on the ionisation data
        return ER*np.ones(len(self.Nuclei())) # all isotopes will have the same ionisation model
    
    def dERdE(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output derivative of ER wrt E_obs, units [keV]/[keV_0]
        """
        deriv = np.gradient(self.ER_samp,self.E_samp) # get the derivative of ER vs E_obs
        return np.interp(E*keV,self.E_samp,deriv)*np.ones(len(self.Nuclei())) # account for both the derivation, and the kg of each isotope per kg of Ge
    
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

        self.si28 = Si28(shell_model)
        self.si29 = Si29(shell_model)
        self.si30 = Si30(shell_model)

        # construct arrays for E_osb --> ER interpolation
        eps = 3.82 #eV
        self.ER_samp = np.arange(0,100*keV,100) # recoil energy in eV range up to 100 keV in steps of 100 eV
        self.E_samp = (Y_Si(self.ER_samp)*volt/eps)*self.ER_samp + self.ER_samp # observed energy in eV

    def Nuclei(self):
        return [[self.si28,0.922],[self.si29,0.047],[self.si30,0.031]]
    
    def ER_E(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output units: [eV] recoil energy keV
        """
        ER = np.interp(E*keV,self.E_samp,self.ER_samp)
        return ER*np.ones(len(self.Nuclei())) # this could alternatively be driven by calibration process
    
    def dERdE(self,E):
        """
        [E] = [keV_0] Observed energy keV

        Output derivative of ER wrt E_obs, units [keV]/[keV_0]
        """
        deriv = np.gradient(self.ER_samp,self.E_samp) # get the derivative of ER vs E_obs
        return np.interp(E*keV,self.E_samp,deriv)*np.ones(len(self.Nuclei()))
    
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
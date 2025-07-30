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

def Y_Ge_Lindhard(ER):
    """
    Ionisation yield for Ge. Taken from the Lindhard assumptions given in arxiv 1304.6773
    [ER] = [eV] recoil energy
    """
    Z = 32
    A = 72.64
    k = 0.133*pow(Z,2/3)*pow(A,-0.5)
    ep = 11.5*ER*pow(Z,-7/3)/keV ## for this we need ER in keV
    g = 3*pow(ep,0.15)+0.7*pow(ep,0.6)+ep
    return k*g/(1+k*g)

def Y_Ge_Sarkis_mean(ER):
    """
    Ionisation yield for Ge using the fit parameters and the corrected Lindhard model
    [ER] = [eV] recoil energy
    """
    k = 0.162
    Z = 32
    A = 72.64
    c0 = 3.0E-4
    c1 = 0.62E-5
    U = 0.02
    cz = 11.5*pow(Z,-7./3.)
    epsR = cz*ER/keV ## for this we need ER in keV
    u = cz*U
    eps = np.where(epsR>= u,epsR - u,0)
    g = 3*pow(eps,0.15) + 0.7*pow(eps,0.6) + eps
    nuL = eps/(1+k*g)
    nu = nuL + c0*pow(eps,0.5) + c1 + u
    return 1-(nu)/(eps+u)

def Y_Ge_Sarkis_min(ER):
    """
    Ionisation yield for Ge using the fit parameters and the corrected Lindhard model
    [ER] = [eV] recoil energy
    """
    Z = 32
    A = 72.64
    k = 0.162 - 0.021
    c0 = (3.0-1.3)*1E-4
    c1 = (0.62-0.12)*1E-5
    U = 0.02-0.015
    cz = 11.5*pow(Z,-7./3.)
    epsR = cz*ER/keV ## for this we need ER in keV
    u = cz*U
    eps = np.where(epsR>= u,epsR - u,0)
    g = 3*pow(eps,0.15) + 0.7*pow(eps,0.6) + eps
    nuL = eps/(1+k*g)
    nu = nuL + c0*pow(eps,0.5) + c1 + u
    return 1-(nu)/(eps+u)

def Y_Ge_Sarkis_max(ER):
    """
    Ionisation yield for Ge using the fit parameters and the corrected Lindhard model
    [ER] = [eV] recoil energy
    """
    Z = 32
    A = 72.64
    k = 0.162 + 0.028
    c0 = (3.0+1.3)*1E-4
    c1 = (0.62+0.12)*1E-5
    U = 0.02+0.01
    cz = 11.5*pow(Z,-7./3.)
    epsR = cz*ER/keV ## for this we need ER in keV
    u = cz*U
    eps = np.where(epsR>= u,epsR - u,0)
    g = 3*pow(eps,0.15) + 0.7*pow(eps,0.6) + eps
    nuL = eps/(1+k*g)
    nu = nuL + c0*pow(eps,0.5) + c1 + u
    return 1-(nu)/(eps+u)

def Y_Ge_LTP(ER):
    """
    Ionisation yield for Ge. Taken from the LTP document
    [ER] = [eV] recoil energy
    """
    ER_samp = [0.0000107,0.000127,0.000153,0.000181,0.000215,0.000256,0.000304,0.000361,0.000428,0.000509,0.000604,0.000718,0.000852,0.00101,0.0012,0.00143,0.0017,0.00201,0.00239,0.00284,0.00335,0.00401,0.00476,0.00565,0.00671,0.00794,0.00948,0.0113,0.0133,0.0155,0.0178,0.0205,0.0227,0.0247,0.027,0.0292,0.0318,0.0349,0.0383,0.0421,0.0462,0.052,0.0603,0.0711,0.0844,0.1,0.119,0.141,0.168,0.199,0.237,0.281,0.334,0.397,0.471,0.56,0.665,0.789,0.937,1.11,1.32,1.57,1.86,2.21,2.63,3.12,3.71,4.41,5.23,6.21,7.38,8.76,10.2,11.9,13.9,16.3,18.9,21.5,24.4,27.9,31.6,35.5,40.2,45.2,50.5,56.3,62.8,69.0,75.2,82.6,90.7,97.3] # keV
    y_samp = [0.00616,0.00614,0.00714,0.00715,0.00715,0.00716,0.00716,0.00716,0.00755,0.00817,0.00818,0.00818,0.00884,0.00919,0.00919,0.0097,0.0102,0.0105,0.0112,0.0119,0.0123,0.0134,0.0147,0.0159,0.0177,0.0197,0.0225,0.026,0.0303,0.0349,0.0409,0.0477,0.0533,0.0587,0.0644,0.0701,0.0765,0.0838,0.091,0.097,0.103,0.109,0.115,0.12,0.124,0.127,0.13,0.133,0.136,0.139,0.142,0.145,0.149,0.152,0.156,0.16,0.163,0.167,0.171,0.175,0.179,0.183,0.187,0.191,0.196,0.2,0.205,0.21,0.215,0.22,0.225,0.231,0.236,0.242,0.247,0.253,0.258,0.264,0.269,0.275,0.28,0.285,0.291,0.297,0.303,0.308,0.314,0.32,0.325,0.331,0.337,0.341] # unitless
    y = np.interp(ER/keV,ER_samp,y_samp) # turn input into keV for interpolation
    return y

def Y_Si_LTP(ER):
    """
    Ionisation yield for Si. Taken from the LTP document
    [ER] = [eV] recoil energy
    """
    ER_samp = [0.000108,0.000128,0.000152,0.00018,0.000214,0.000254,0.000302,0.000359,0.000426,0.000506,0.000601,0.000714,0.000848,0.00101,0.0012,0.00142,0.00169,0.002,0.00238,0.00283,0.00336,0.00398,0.00473,0.00562,0.00668,0.00793,0.00942,0.0112,0.0133,0.0158,0.0187,0.0223,0.0264,0.0314,0.0373,0.0443,0.0526,0.0624,0.0742,0.0881,0.105,0.124,0.148,0.175,0.208,0.247,0.294,0.349,0.414,0.492,0.571,0.652,0.75,0.87,1.03,1.21,1.4,1.61,1.84,2.09,2.37,2.7,3.09,3.5,3.97,4.49,5.09,5.73,6.39,7.13,7.89,8.73,9.75,10.9,12.1,13.6,15.6,17.9,20.6,23.6,26.7,30.3,34.3,39.2,43.7,48.4,53.2,58.4,64.1,70.4,76.8,83.0,89.7,96.3] # keV
    y_samp = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000997,0.00242,0.00438,0.00626,0.00876,0.011,0.0135,0.0166,0.0199,0.024,0.0284,0.0335,0.0386,0.0449,0.0519,0.06,0.0693,0.0778,0.0871,0.0962,0.106,0.115,0.124,0.133,0.142,0.152,0.162,0.172,0.182,0.191,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.269,0.279,0.29,0.3,0.311,0.32,0.331,0.34,0.35,0.36,0.371,0.381,0.392,0.404,0.415,0.425,0.435,0.445,0.456,0.466,0.477,0.486,0.496,0.504] # unitless
    y = np.interp(ER/keV,ER_samp,y_samp) # turn input into keV for interpolation
    return y


class GeHV(Detector):
    def __init__(self, volt, shell_model="Fitz",y_model=Y_Ge_LTP,A=None,B=None,sig_E=None):
        """
        SuperCDMS-like Ge HV detector. 
        Input parameters:
        -----------------
        volt = voltage detector is held at [V]
        shell_model = nuclear shell model assumptions; Fitz, JUN45, jj44b
        y_model = ionisation yield function. Can use one of the ones defined above, or pass any function that takes a single argument ER
        A = resolution model parameter that scales with E^2 [unitless]
        B = resolution model parameter that scales with E [keV]
        sig_E = resolution model parameter independent of E [keV]
        """
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
        self.ER_samp = np.arange(0,50*keV,100) # recoil energy in eV range up to 100 keV in steps of 100 eV
        self.E_samp = (y_model(self.ER_samp)*volt+eps)*self.ER_samp/(eps+volt)# observed energy in eV

        # Parameters for resolution
        if A==None:
            self.A = 5E-3
        else:
            self.A = A
        if B==None:
            self.B = 0.7/keV # 0.7 eV
        else:
            self.B = B
        if sig_E==None:
            self.sig_E = 10/keV # 10 eV
        else:
            self.sig_E = sig_E

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
        return [0.01,30]
    
    def Emax(self):
        return 20
    
    def DeltaE(self,E):
        # Lets assume CDMSlite values: arxiv 1911.11905. These are for a Ge iZIP, so we'll at least adjust the baseline res to what we hope for
        return np.sqrt(self.B*E+pow(self.A*E,2)+pow(self.sig_E,2)) # these resolutions are defined for eV so need to convert to keV for comp of DeltaE
    
    def Res(self,E1,E2):
        # We assume E1 is the observed energy (E' in accompanying documentation) and E2 is the energy that will be integrated over (E_ee in accompanying documentation)
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        return 0.85
    
class SiHV(Detector):
    def __init__(self, volt, shell_model="Fitz", y_model=Y_Si_LTP,A=None,B=None,sig_E=None):
        """
        SuperCDMS-like Si HV detector. 
        Input parameters:
        -----------------
        volt = voltage detector is held at [V]
        shell_model = nuclear shell model assumptions; Fitz, USDB
        y_model = ionisation yield function. Can use one of the ones defined above, or pass any function that takes a single argument ER
        A = resolution model parameter that scales with E^2 [unitless]
        B = resolution model parameter that scales with E [keV]
        sig_E = resolution model parameter independent of E [keV]
        """
        ## allow for initialisation with different shell models
        self.shell_model = shell_model
        self.V = volt # voltage detector is run at in volts

        self.si28 = Si28(shell_model)
        self.si29 = Si29(shell_model)
        self.si30 = Si30(shell_model)

        # construct arrays for E_osb --> ER interpolation
        eps = 3.82 #eV
        self.ER_samp = np.arange(0,100*keV,100) # recoil energy in eV range up to 100 keV in steps of 100 eV
        self.E_samp = (y_model(self.ER_samp)*volt+eps)*self.ER_samp/(eps+volt)# observed energy in eV

        # Parameters for resolution
        if A==None:
            self.A = 5E-3
        else:
            self.A = A
        if B==None:
            self.B = 0.7/keV # 0.7 eV
        else:
            self.B = B
        if sig_E==None:
            self.sig_E = 13/keV # 34 eV
        else:
            self.sig_E = sig_E

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
        return [0.1,10]
    
    def Emax(self):
        return 20
    
    def DeltaE(self,E):
        # Lets assume CDMSlite values: arxiv 1911.11905. These are for a Ge iZIP, so we'll at least adjust the baseline res to what we hope for
        return np.sqrt(self.B*E+pow(self.A*E,2)+pow(self.sig_E,2)) # these resolutions are defined for eV so need to convert to keV for comp of DeltaE
    
    def Res(self,E1,E2):
        # We assume E1 is the observed energy (E' in accompanying documentation) and E2 is the energy that will be integrated over (E_ee in accompanying documentation)
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        return A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
    
    def Eff(self,E):
        return 0.85
from neutrinomodel import NeutrinoModel
import numpy as np
from constants import *
import scipy.integrate as integrate
import scipy

######################################################
### Class definition for standard Supernova neutrino background model
#### https://doi.org/10.1016/j.astropartphys.2023.102890

class SupernovaNeutrino(NeutrinoModel):
    def Flux_neutrino(self,Ev,source,epsilon=1,E_ave=12*MeV,gamma=3,L=1):
        """
        Calculate neutrino flux at the detector from the supernova.
    
        Parameters:
        - Ev: Energy of neutrinos [eV]
	- flavor: e for electron neutrino, ea for electron anti-neutrino, x for the rest
        - L: Distance from the supernova to the detector (Earth) [cm]
        - Ni: The number of emitted neutrinos for flavor i, data from https://www.sciencedirect.com/science/article/pii/S0927650512001132
        - Ti: The temperature of emitted neutrinos for flavor i
    
        Returns:
        - phi(E) in [cm-2eV-1]
        """
        if source == 'supernova_e' or source == 'supernova_ea'or source == 'supernova_mu' or source == 'supernova_mua' or source == 'supernova_tau' or source == 'supernova_taua' :
            D=L*10*kpc_to_cm
            epsilon=epsilon*5e52*erg
            result = 1/(4*np.pi*D**2)*epsilon/E_ave*Ev**gamma/scipy.special.gamma(1+gamma)*((1+gamma)/E_ave)**(1+gamma)*np.exp(-(1+gamma)*Ev/E_ave)
        else:
            print("Please choose the right flavor of neutrinos. Options are supernova_e,supernova_ea,supernova_mu,supernova_mua,supernova_tau,supernova_taua")
            result=0
        return result
    
            
    def dRdER(self, Ed, channel=None, source=None,
              Ev=np.logspace(-2,3,501)*MeV, E_ave=12*MeV):
        """
            [Ev]      = [eV]  Neutrino energy
            [Ed]      = [eV]  Detected energy
            [E_ave]   = [eV]  Average energy of Neutrino
            [channel] = [1]   Neutrino interaction channel in the detector
            [source]  = []    species of the neutrino
        """
        if channel == "IBD": 
            ##inverse beta decay cross section (IBD) https://arxiv.org/pdf/astro-ph/0302055.pdf
            def Sig(Ev, Eout, channel,source):
                pep=np.sqrt(Eout**2-me**2)
                Ev=Ev/MeV
                Eout=Eout/MeV
                pep=pep/MeV
                cross_sec = 1e-43*Eout*pep*Ev**(-0.07056+0.02018*np.log(Ev)-0.001953*(np.log(Ev))**3)
                return cross_sec
            Eout=Ed-me
            Ev=Eout+Delta_np
            if Ev > Delta_np+me:
                integral = Sig(Ev,Eout,channel,source) * self.Flux_neutrino(Ev,source,1,E_ave,3,1) * MeV
            else:
                integral = 0
                
                
        elif channel == "pES":
            def Sig(Ev,Eout,channel,source):
                cross_sec = 4.83e-42/MeV*(1+466*Eout/MeV/(Ev/MeV)**2)
                return cross_sec
            def integrand(Ev,Ed,e_ave,channel):
                return Sig( Ev, Ed, channel,source) * self.Flux_neutrino(Ev,source,1,e_ave,3,1)
            Emin=np.sqrt(Ed*mproton/2)
            integral = integrate.quad(lambda Evl:integrand(Evl,Ed,E_ave,channel),Emin,max(Ev))[0] * MeV

            
        elif channel == "eES":
            def Sig(Ev, Eout, channel,source):
                if source == 'supernova_e' or source == 'supernova_ea':
                    epsilonN=-0.5-sin2thetaW
                    epsilonP=-sin2thetaW
                else:
                    epsilonN=0.5-sin2thetaW
                    epsilonP=-sin2thetaW
                if source == 'supernova_e' or source == 'supernova_mu' or source == 'supernova_tau':
                    cross_sec = 2*me*Gf**2/np.pi*(epsilonN**2+epsilonP**2*(1-Eout/Ev)**2-epsilonN*epsilonP*me*Eout/Ev**2)*0.3894*1e-45
                else:
                    cross_sec = 2*me*Gf**2/np.pi*(epsilonP**2+epsilonN**2*(1-Eout/Ev)**2-epsilonN*epsilonP*me*Eout/Ev**2)*0.3894*1e-45
                return cross_sec
            def integrand(Ev,Ed,e_ave,channel,source):
                return Sig( Ev, Ed, channel,source) * self.Flux_neutrino(Ev,source,1,e_ave,3,1)
            Emin=Ed/2+np.sqrt(Ed*(Ed+2*me))/2
            integral = integrate.quad(lambda Evl:integrand(Evl,Ed,E_ave,channel,source),Emin,max(Ev))[0] * MeV

            
        elif channel == 'C12_NC':##https://www.sciencedirect.com/science/article/pii/0370157379900103
            def Sig(Ev,channel,source):
                sig=1.08e-38*((Ev - w_value_C12_M1) / M_nucleus)**2
                return sig
            def integrand(Ev,e_ave,channel):
                return Sig(Ev,channel,source)* self.Flux_neutrino(Ev,source, 1,e_ave,3,1)
            Emin = w_value_C12_M1
            integral = integrate.quad(lambda Evl:integrand(Evl,E_ave,channel),Emin,max(Ev))[0] * branching_ratio_C12_M1

        return integral

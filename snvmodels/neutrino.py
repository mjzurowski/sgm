from neutrinomodel import NeutrinoModel
import numpy as np
from constants import *
from targets.c12 import C12
import scipy.integrate as integrate
from snvmodels.solarndists import SolarNDist
import sympy
import scipy 

x = sympy.symbols('x')

######################################################
### Class definition for standard Supernova neutrino background model
#### https://doi.org/10.1016/j.astropartphys.2023.102890

class Neutrino(NeutrinoModel):
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
#        elif source == 'supernova_total':
#            result = 
        elif source == "solar_f17":
            f17NDist = SolarNDist("f17",0.0005912e10)
            result=f17NDist.dist(Ev/MeV)/MeV*24*3600
        elif source == "solar_o15":
            o15NDist = SolarNDist("o15",0.05031e10)
            result=o15NDist.dist(Ev/MeV)/MeV*24*3600 
        elif source == "solar_n13":
            n13NDist = SolarNDist("n13",0.05712e10)
            result=n13NDist.dist(Ev/MeV)/MeV*24*3600
        elif source == "solar_hep":
            hepNDist = SolarNDist("hep",0.000000788e10)
            result=hepNDist.dist(Ev/MeV)/MeV*24*3600
        elif source == "solar_b8":
            b8NDist = SolarNDist("b8",0.0005822e10)
            result=b8NDist.dist(Ev/MeV)/MeV*24*3600
        elif source == "solar_pp":
            ppNDist = SolarNDist("pp",5.938e10)
            result=ppNDist.dist(Ev/MeV)/MeV*24*3600
        elif source == "solar_be7":
            result = sympy.DiracDelta(x-384.3e3)* 0.103*0.4857e10 + sympy.DiracDelta(x-861.8e3)* 0.897*0.4857e10
        elif source =='solar_pep':
            result = sympy.DiracDelta(x-1442e3)*0.01401e10
        else:
            print("Please choose the right flavor of neutrinos. Options are supernova_e,supernova_ea,supernova_mu,supernova_mua,supernova_tau,supernova_taua")
            result=0
        return result   
    
            
    def dRdEd(self,Ev,Ed,E_ave,channel,source):
        """
            [Ev]      = [eV]  Neutrino energy
            [Ed]      = [eV]  Detected energy
            [E_ave]   = [eV]  Average energy of Neutrino
            [channel] = [1]   Neutrino interaction channel in the detector
            [source]  = []    species of the neutrino
        """
        # Intergration over neutrino spectrum
        #proton, C12 and electron number density in LAB
        Np=7e28 #per ton from https://arxiv.org/pdf/1712.06985.pdf
        NC=4.4e28
        Ne=Np+6*NC
        
        if channel == "IBD":
            ##inverse beta decay cross section (IBD) https://arxiv.org/pdf/astro-ph/0302055.pdf
            def Sig(Ev, Eout, channel,source):
                pep=np.sqrt(Eout**2-me**2)
                cross_sec = 9.52e-44*(Eout*pep/MeV**2)
                return cross_sec
            Eout=Ed-me
            Ev=Eout+Delta_np
            if Ev > Delta_np+me:
                integral = Np * Sig(Ev,Eout,channel,source) * self.Flux_neutrino(Ev,source,1,E_ave,3,1)
            else:
                integral = 0
                
        elif channel == "pES":
            def Sig(Ev,Eout,channel,source):
                cross_sec = 4.83e-42/MeV*(1+466*Eout/MeV/(Ev/MeV)**2)
                return cross_sec
            def integrand(Ev,Ed,e_ave,channel):
                return Np * Sig( Ev, Ed, channel,source) * self.Flux_neutrino(Ev,source,1,e_ave,3,1)
            Emin=np.sqrt(Ed*mproton/2)
            integral = integrate.quad(lambda Evl:integrand(Evl,Ed,E_ave,channel),Emin,max(Ev))[0]
            
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
            def integrand(Ev,Ed,e_ave,channel):
                return Ne * Sig( Ev, Ed, channel,source) * self.Flux_neutrino(Ev,source,1,e_ave,3,1)
            Emin=Ed/2+np.sqrt(Ed*(Ed+2*me))/2
            integral = integrate.quad(lambda Evl:integrand(Evl,Ed,E_ave,channel),Emin,max(Ev))[0]
            
        elif channel == 'C12_NC':
            w=15.1*MeV
            beta=1
            kai=1
            def Sig(Ev,channel,source):
                sig=1.08e-38*((Ev-w)/M_nucleus)**2*beta**2*kai**2
                return sig
            def integrand(Ev,e_ave,channel):
                return NC * Sig(Ev,channel,source)* self.Flux_neutrino(Ev,source, 1,e_ave,3,1)
            Emin=w
            integral = integrate.quad(lambda Evl:integrand(Evl,E_ave,channel),Emin,max(Ev))[0]

        elif channel == 'C12_CEvNS':
            def Sig(Target,Ev,ER):
                FF = Target.Helm(ER)**2 ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
                cross_sec = Gf**2*Target.mT()*1e-9/(8*np.pi)*(Target.Z()*(4*sin2thetaW-1)+(Target.A()-Target.Z()))**2*(2-ER*Target.mT()/Ev**2)*FF*0.389379*1e-27/1e9 #cm2/eV
                return cross_sec
            def integrand(Ev,ER,Target,source):
                return 4.4e28 * Sig(Target, Ev, ER) * self.Flux_neutrino(Ev,source)
            Target=C12()
            integral = integrate.quad(lambda Evl:integrand(Evl,Ed,Target,source),np.sqrt(Ed*Target.mT()/2),max(Ev))[0]
            
        return integral



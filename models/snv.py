from numodel import NuModel
import numpy as np
from constants import *
import scipy.integrate as integrate

######################################################
### Class definition for standard Supernova neutrino background model
#### https://doi.org/10.1016/j.astropartphys.2023.102890

class SNv(NuModel):
    def Flux_neutrino(self,Ev,source):
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
        L = 3.086e+22 #10 kpc [cm]
        Ni_e = 3e57
        T_e = 3.5*MeV
        Ni_ea = 2.1e57
        T_ea = 5*MeV
        Ni_x = 5.2e57
        T_x = 8*MeV 
        if source == 'supernova_e':
            result = Ni_e*Ev**2/(4*np.pi*L**2*2*T_e**3)*np.exp(-Ev/T_e)
        elif source == 'supernova_ea':
            result = Ni_ea*Ev**2/(4*np.pi*L**2*2*T_ea**3)*np.exp(-Ev/T_ea)
        elif source == 'supernova_x':
            result = Ni_x*Ev**2/(4*np.pi*L**2*2*T_x**3)*np.exp(-Ev/T_x)
        elif source == 'supernova':
            result = Ni_e*Ev**2/(4*np.pi*L**2*2*T_e**3)*np.exp(-Ev/T_e) +Ni_ea*Ev**2/(4*np.pi*L**2*2*T_ea**3)*np.exp(-Ev/T_ea)+Ni_x*Ev**2/(4*np.pi*L**2*2*T_x**3)*np.exp(-Ev/T_x)
        else:
            print("Please choose the right flavor of neutrinos. Options are 'e' for electron neutrino, 'ea' for electron antineutrino and 'x' for the rest.")
            result=0
        return result   
    
    def Emin(self,Target,ER):
       return np.sqrt(ER*Target.mT()/2)

    def Sig(self, Target, Ev, ER):
        Gf = 1.166364*1e-5 ##Fermi constant in natural unit GeV
        sin2thetaW = 0.2229 ## Weak mixing angle from https://physics.nist.gov/cgi-bin/cuu/Value?sin2th
        FF = Target.Helm(ER)**2 ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = Gf**2*Target.mT()*1e-9/(8*np.pi)*(Target.Z()*(4*sin2thetaW-1)+(Target.A()-Target.Z()))**2*(2-ER*Target.mT()/Ev**2)*FF*0.389379*1e-27/1e9 #cm2/eV
        return cross_sec
            
    
    def dRdER(self, Target, ER, Ev, source):
        # Intergration over neutrino spectrum
        def integrand(Ev,ER,Target,source):
            return Target.N_T() * self.Sig(Target, Ev, ER) * self.Flux_neutrino(Ev,source)
        #integral = np.array([integrate.quad(lambda Evl:integrand(Evl,ERi),self.Emin(Target,ERi),max(Ev),limit=int(1E8))[0] for ERi in ER]) ## this integral could probs be optimised
        if source in ['supernova_e','supernova_ea','supernova_x','supernova']:
            integral = integrate.quad(lambda Evl:integrand(Evl,ER,Target,source),self.Emin(Target,ER),max(Ev))[0]  ## this integral could probs be optimised
        
        #integral = np.trapz(integrand(Ev,ER)[Ev>self.Emin(Target,ER)],Ev[Ev>self.Emin(Target,ER)])
        #print(integral)
        return integral * keV

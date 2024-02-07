from snvmodel import SNvModel
import numpy as np
from constants import *
import scipy.integrate as integrate

######################################################
### Class definition for standard Supernova neutrino background model
#### https://doi.org/10.1016/j.astropartphys.2023.102890

class SNv(SNvModel):
    def Flux_neutrino(self,Ev,flavor):
        """
        Calculate neutrino flux at the detector from the supernova.
    
        Parameters:
        - Ev: Energy of neutrinos
        - L: Distance from the supernova to the detector (Earth)
        - Ni: The number of emitted neutrinos for flavor i, data from https://www.sciencedirect.com/science/article/pii/S0927650512001132
        - Ti: The temperature of emitted neutrinos for flavor i
    
        Returns:
        - phi(E)
        """
        L = 3.086e+22 #10 kpc [cm]
        if flavor == 'e':
            Ni = 3e57
            T = 3.5e3
        elif flavor == 'ea':
            Ni = 2.1e57
            T = 5e3 
        elif flavor == 'x':
            Ni = 5.2e57
            T = 8e3 
        else:
            print("Please choose the right flavor of neutrinos. Options are 'e' for electron neutrino, 'ea' for electron antineutrino and 'x' for the rest.")
        
        result = Ni*Ev**2/(4*np.pi*L**2*2*T**3)*np.exp(-Ev/T)
        return result   
    
    def Emin(self,Target,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return np.sqrt(ER*Target.mT()*1e-3/2)

    def Sig(self,Target,Ev,ER):
        """
        Differential cross section as a function of recoil energy ER in [KeV]
        Target: target nucleus
        mX: DM mass [eV]
        ER: recoil energy [eV]
        """
        Gf = 1.166364*1e-5 ##Fermi constant in natural unit GeV
        sin2thetaW = 0.2229 ## Weak mixing angle from https://physics.nist.gov/cgi-bin/cuu/Value?sin2th
        FF = Target.Helm(ER*1e3)**2 ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = Gf**2*Target.mT()*1e-9/(8*np.pi)*(Target.Z()*(4*sin2thetaW-1)+(Target.A()-Target.Z()))**2*(2-ER*Target.mT()*1e-3/Ev**2)*FF*0.389379*1e-27/1e6 ##cm2/keV
        return cross_sec
            
    
    def dRdER(self,Target,Ev,ER):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        Flux_Total=self.Flux_neutrino(Ev,'e')+self.Flux_neutrino(Ev,'ea')+self.Flux_neutrino(Ev,'x')
        #print('cpd: ', cpd_conversion)
        #print('Target: ',Target.N_T())
        #print('Sig: ', self.Sig(Target,Ev,ER))
        #print('Flux: ', Flux_Total)
        # Intergration over neutrino spectrum
        def integrand(Ev,ER):
            return Target.N_T() * self.Sig(Target, Ev, ER) * (self.Flux_neutrino(Ev,'e')+self.Flux_neutrino(Ev,'ea')+self.Flux_neutrino(Ev,'x'))
        #integral = np.array([integrate.quad(lambda Evl:integrand(Evl,ERi),self.Emin(Target,ERi),max(Ev),limit=int(1E8))[0] for ERi in ER]) ## this integral could probs be optimised
        integral = integrate.quad(lambda Evl:integrand(Evl,ER),self.Emin(Target,ER),max(Ev))[0]  ## this integral could probs be optimised
        #integral = np.trapz(integrand(Ev,ER)[Ev>self.Emin(Target,ER)],Ev[Ev>self.Emin(Target,ER)])
        #print(integral)
        return integral

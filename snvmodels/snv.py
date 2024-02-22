from snvmodel import SNvModel
import numpy as np
from constants import *
import scipy.integrate as integrate
from snvmodels.solarndists import SolarNDist
import sympy

x = sympy.symbols('x')

######################################################
### Class definition for standard Supernova neutrino background model
#### https://doi.org/10.1016/j.astropartphys.2023.102890

class SNv(SNvModel):
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
            print("Please choose the right flavor of neutrinos. Options are 'e' for electron neutrino, 'ea' for electron antineutrino and 'x' for the rest.")
            result=0
        return result   
    
    def Emin(self,Target,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [eV]
       """
       return np.sqrt(ER*Target.mT()/2)

    def Sig(self,Target,Ev,ER):
        """
        Differential cross section as a function of recoil energy ER in [eV]
        Target: target nucleus
        mX: DM mass [eV]
        ER: recoil energy [eV]
        Sig: cross section [cm2/eV]
        """
        Gf = 1.166364*1e-5 ##Fermi constant in natural unit GeV
        sin2thetaW = 0.2229 ## Weak mixing angle from https://physics.nist.gov/cgi-bin/cuu/Value?sin2th
        FF = Target.Helm(ER)**2 ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = Gf**2*Target.mT()*1e-9/(8*np.pi)*(Target.Z()*(4*sin2thetaW-1)+(Target.A()-Target.Z()))**2*(2-ER*Target.mT()/Ev**2)*FF*0.389379*1e-27/1e9 #cm2/eV
        return cross_sec
            
    
    def dRdER(self,Target,Ev,ER,source):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm2/eV] cross section

        Output units: [cm^2]/[eV] 
        """
        # Intergration over neutrino spectrum
        def integrand(Ev,ER,Target,source):
            return Target.N_T() * self.Sig(Target, Ev, ER) * self.Flux_neutrino(Ev,source)
        #integral = np.array([integrate.quad(lambda Evl:integrand(Evl,ERi),self.Emin(Target,ERi),max(Ev),limit=int(1E8))[0] for ERi in ER]) ## this integral could probs be optimised
        if source in ['supernova_e','supernova_ea','supernova_x','supernova']:
            integral = integrate.quad(lambda Evl:integrand(Evl,ER,Target,source),self.Emin(Target,ER),max(Ev))[0]  ## this integral could probs be optimised
        elif source in ['solar_be7','solar_pep']:
            integral = sympy.integrate(integrand(x,ER,Target,source),(x,self.Emin(Target,ER),float('inf')))
        elif source in ['solar_f17','solar_o15','solar_n13','solar_hep','solar_b8','solar_pp']:
            integral = np.trapz(integrand(Ev[Ev>self.Emin(Target,ER)],ER,Target,source),Ev[Ev>self.Emin(Target,ER)])
        
        #integral = np.trapz(integrand(Ev,ER)[Ev>self.Emin(Target,ER)],Ev[Ev>self.Emin(Target,ER)])
        #print(integral)
        return integral

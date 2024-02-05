from dmmodel import DMModel
import numpy as np
from constants import *

### Class definition for MIGDAL Effect using Standard WIMP Assumptions.

class MIGDAL(DMModel):

    def vmin(self, Target, mX, ER, EE):
        
        """
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [EE] = [eV] Ionization energy (Electron + Binding)

        Output units: [km/s]

        This should return an array of vmin values accounting for each of the transition elements.x
        """
        
        vmins = [] # VMIN is returned as a list based on the transition energies in the target.
        for Eb in Target.eNL():
            vmins.append(c*1e-3*(Target.mT()*ER+Target.mu_T(mX)*(-1*Eb+EE))/(Target.mu_T(mX)*np.power(2*Target.mT()*ER,0.5)))
        
        return vmins
    
 
    def dRdER(self, Target, mX, ER, sig, VelDist, EE):
        
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX]  = [eV] DM mass
        [ER]  = [eV] DM recoil energy
        [EE]  = [eV] Ionization energy (Electron + Binding)
        [sig] = [cm]^2 cross section
        [EE]  = [eV] Ionization energy (Electron + Binding)

        Output units: [cm^2]/[eV] 
        """
        FF = Target.Helm(ER)
        PreFactor = (sig)*(VelDist.rho)*(Target.mT()**2)*FF/(2*mX*Target.mu_N(mX)**2)
        
        # vmins  = [vmin for vmin in self.vmin(Target, ER, mX, EE)] # These will give an array of vmins for the given EE and Transition energies in the target.
        vmins = self.vmin(Target, mX, ER, EE)        
        gdists = VelDist.gdist(vmins) # Similarly taking and evaluating for the relevant vmin.
        probs  = Target.eTrans(EE)

        # print('vmin', vmins)
        # print('gdist', gdists)
        # print('prob', probs)


        summed_contribution = np.sum(gdists*probs)
        
        # print(summed_contribution)



        return cpd_conversion*PreFactor*np.array(summed_contribution)


""" Requre the addition of both a probability calculation and ionization rate calculation."""
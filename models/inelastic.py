from dmmodel import DMModel
import numpy as np
from constants import *

### Class definition for inelastic SI WIMP.

class SIInel(DMModel):
    def vmin(self,Target,mX,ER,delta):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return (c*1.E-3)*np.abs((Target.mT()*ER/Target.mu_T(mX))+delta)/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist,delta):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        FF = Target.Helm(ER)**2 ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        dsigdER = (1/kg_to_eV)*sig*FF*Target.A()*Target.A()/(2*Target.N_T()*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER,delta)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)

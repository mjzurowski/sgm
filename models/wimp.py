from dmmodel import DMModel
import numpy as np
from constants import *

### Class definition for standard SI WIMP.

class SIWIMP(DMModel):
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return (c*1.E-3)*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dsigdER(self,Target,mX,ER,sig):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        FF = Target.F11(ER,1/np.sqrt(2),1/np.sqrt(2)) ## cross section with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = sig*2 ## as we normalise the nucleon vector [cn, cp] our built in cross section is sigma_N. Multiply by 2 to get this from sigma_p (which is generally assumed to be the input)
        return cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX))

######################################################
### Class definition for standard SI WIMP with Helm form factors
    
class SIWIMP_Helm(DMModel):
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return (c*1.E-3)*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dsigdER(self,Target,mX,ER,sig):
        """
        For this model, we take the Helm Form Factor rather than those defined from nrefts. Note that it should be ~ equal to the O1 SIWIMP
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        FF = Target.Helm(ER)**2 ## cross section with couplings. Note that proton and neutron couplings are normalised to 1
        return (1/kg_to_eV)*sig*FF*Target.A()*Target.A()/(2*Target.N_T()*Target.mu_N(mX)*Target.mu_N(mX))
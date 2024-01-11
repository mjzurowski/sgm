from dmmodel import DMModel
import numpy as np

c = 299792458. # speed of light in m/s

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
        return sig*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX))

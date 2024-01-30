from dmmodel import DMModel
import numpy as np
from constants import *

### Class definition for MIGDAL Effect using Standard WIMP Assumptions.

class MIGDAL(DMModel):
    def vmin(self,Target,mX,ER):
        """
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: [km/s]
        """
        return (c*1.E-3)*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def transitionprob():
        pass
    
    def dsigdER(self, Target, mX, ER, sig):
        pass


""" Requre the addition of both a probability calculation and ionization rate calculation."""
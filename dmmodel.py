from abc import ABC, abstractmethod
from constants import *
import numpy as np

"""
This abstract class can be used to define general (elastic) interaction models constructed from some differential cross section
Note that this is not the ONLY way to define models, but is most frequently used and so most convenient to make an abtract class
Examples for DM models not based on this class can be found in models/
"""

class DMModel(ABC):
    def vmin(self,Target,mX,ER,**kwargs):
        """
        Minimum velocity that can produce recoil energy ER in [km/s]
        [Target]: target nucleus
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: [km/s]
        """
        return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)

    # @abstractmethod
    # def dsigdER(self,Target,ER,mX,sig):
    #     """
    #     Differential cross section in units of [cm^2]/[eV]
    #     Note the exact form of this will depend on the assumptions being used for the form factors
    #     """
    #     pass

    def dRdER(self,Target,ER,mX,sig,VelDist,**kwargs):
        """
        Interaction rate as a function of recoil energy in counts/[day]/[kg]/[keV]
        Inputs:
            Target: target nucleus (see target.py)
            mX: DM mass [eV]
            ER: recoil energy [eV]
            sig: DM cross section [cm]^2
            dist: velocity distribution [unitless]
        """
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*VelDist.gdist(vm)*self.dsigdER(Target,ER,mX,sig)



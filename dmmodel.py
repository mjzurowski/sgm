from abc import ABC, abstractmethod
from constants import *

class DMModel(ABC):
    @abstractmethod
    def vmin(self,Target,mX,ER):
        """
        Minimum velocity that can produce recoil energy ER in [km/s]
        Target: target nucleus
        mX: DM mass [eV]
        ER: recoil energy [eV]
        """
        pass

    @abstractmethod
    def dRdER(self,Target,mX,ER,sig,VelDist):
        """
        Interaction rate as a function of recoil energy in counts/[day]/[kg]/[keV]
        Inputs:
            Target: target nucleus (see target.py)
            mX: DM mass [eV]
            ER: recoil energy [eV]
            sig: DM cross section [cm]^2
            dist: velocity distribution [unitless]
        """
        pass



from abc import ABC, abstractmethod
from constants import *

class NuModel(ABC):
    @abstractmethod
    def Flux_neutrino(self, Ev, flavor):
        """
        Calculate neutrino flux at the detector from the supernova.
    
        Parameters:
        - Ev: Energy of neutrinos with a certain flavor.
        - flavor: The flavor of neutrinos
        Returns:
        - phi(E)
        """
        pass
        
    @abstractmethod    
    def Emin(self,Target,mX,ER):
        """
        Minimum Energy that can produce recoil energy ER in [MeV]
        Target: target nucleus
        mX: DM mass [eV]
        ER: recoil energy [eV]
        """
        pass

    @abstractmethod    
    def Sig(self,Target,mX,Ev, ER):
        """
        Differential cross section as a function of recoil energy ER in [KeV]
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



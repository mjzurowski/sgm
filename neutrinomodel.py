from abc import ABC, abstractmethod
from constants import *

class NeutrinoModel(ABC):
    @abstractmethod
    def Flux_neutrino(self,Ev,source,epsilon=1,E_ave=12*MeV,gamma=3,L=1):
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
    def dRdEd(self,Ev,Ed,E_ave,channel,source):
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


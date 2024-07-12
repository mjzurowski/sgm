from abc import ABC, abstractmethod
import numpy as np
from constants import *

#### Fundamental properties of a target that won't change between detectors (e.g., mass, form factors)
##### Note that a lot of these functions depend on additional data. These are included in subfolders in targets/, where info.txt files can be found that give detail on their units and source

class Target(ABC):
    def spin_dep(self, jx):
        return (4*jx/3)*(jx+1) #### double check this defn as I think it should depend on the target
    
    @abstractmethod
    def A(self):
        """
        Mass number of target [unitless]
        """
        pass

    def N_T(self):
        """
        Target density [atoms/kg]
        """
        return 6.02E26/self.A()

    def mT(self):
        """
        Mass of target [eV]
        """
        return self.A()*mp
    
    def mu_T(self,mX):
        """
        Reduced mass of system with DM [eV]
        """
        return self.mT()*mX/(self.mT()+mX)
    
    def mu_N(self,mX):
        """
        Reduced mass of nucleon with DM [eV]
        """
        return mp*mX/(mp+mX)

    def B(self):
        """
        Units of [eV]^-1
        """
        return np.power(41.467/(45.*np.power(self.A(),-1./3.)-25.*np.power(self.A(),-2./3.)),0.5)
    
    def Q(self,ER):
        """
        Units of [eV]
        """
        return np.power(2.*self.mT()*ER, 0.5)
    
    def Y(self,ER):
        """
        Unitless parameter used to define form factors
        """
        return np.power((1/(197.327 *1E6))*self.Q(ER)*self.B()/2.,2.)
    
    def FMpp(self,ER):
        pass

    def FMnp(self,ER):
        pass

    def FMpn(self,ER):
        pass

    def FMnn(self,ER):
        pass

    def FS1pp(self,ER):
        pass

    def FS1np(self,ER):
        pass

    def FS1pn(self,ER):
        pass

    def FS1nn(self,ER):
        pass

    def FS2pp(self,ER):
        pass

    def FS2np(self,ER):
        pass

    def FS2pn(self,ER):
        pass

    def FS2nn(self,ER):
        pass
        
    def FPhi2pp(self,ER):
        pass

    def FPhi2np(self,ER):
        pass

    def FPhi2pn(self,ER):
        pass

    def FPhi2nn(self,ER):
        pass

    def FPhipp(self,ER):
        pass

    def FPhinp(self,ER):
        pass

    def FPhipn(self,ER):
        pass

    def FPhinn(self,ER):
        pass    
    
    def FDpp(self,ER):
        pass

    def FDnp(self,ER):
        pass

    def FDpn(self,ER):
        pass

    def FDnn(self,ER):
        pass

    def FMPhi2pp(self,ER):
        pass

    def FMPhi2np(self,ER):
        pass

    def FMPhi2pn(self,ER):
        pass

    def FMPhi2nn(self,ER):
        pass
    
    def FS1Dpp(self,ER):
        pass

    def FS1Dnp(self,ER):
        pass

    def FS1Dpn(self,ER):
        pass

    def FS1Dnn(self,ER):
        pass

    def Helm(self,ER):
        c1 = 1.23*np.power(self.A(),1/3)-0.6
        s = 0.9
        rn = np.power(c1*c1+(7/3)*0.52*0.52*np.pi*np.pi-5*s*s,0.5) #fm
        q = self.Q(ER)/(197.327*1E6) # Q in units of inverse fm
        return 3*np.exp(-0.5*(s*q)**2)*(np.sin(q*rn)-q*rn*np.cos(q*rn))/(np.power(q*rn,3))
    
##### Electron transition probs
    def eNL(self):
        """
        List of transition energies for the target
        Should have units of eV
        """
        pass

    def eTransE_E(self):
        """
        Files that give electron transition probability as a function of the kinetic energy of the emitted electron E_e (E_e = E_EM - nuclear bonding)
        Note that targets will have multiple of these depending on the transition energies allowed
        These should have units of [eV]^-1 and be passed as a list the same length as eNL
        """
        pass

    def eTrans(self,E_EM):
        """
        Electron transition probability as a function of total electronic energy seen in the detector (E_EM = electron kinetic energy + nuclear bonding). Transforms eTransE_E based on the transition energies
        Should still be a list of the same length (and in units of [eV]^-1) as it later needs to be paired with appropriate vmins
        """
        trans_list = []
        for i in range(0,len(self.eNL)):
            trans_list.append(np.interp(E_EM-self.eNL[i],self.eTransE_E[i][:,0],self.eTransE_E[i][:,1])*np.heaviside(E_EM-self.eNL[i],1))
        return trans_list
    
##### Photoelectric absorption
    def sigma_PE(self,E_gam):
        """
        Photoelectric absorption as a function of photon energy (E_gam) [eV]
        Output units are [cm2/atom]
        """
        pass

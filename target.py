from abc import ABC, abstractmethod
import numpy as np

mp =  0.9314941*1E9 # mass of a nulceon in eV.
#### Fundamental properties of a target that won't change between detectors (e.g., mass, form factors)

class Target(ABC):
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
    
    @abstractmethod
    def FMpp(self,ER):
        pass

    @abstractmethod
    def FMnp(self,ER):
        pass

    @abstractmethod
    def FMpn(self,ER):
        pass

    @abstractmethod
    def FMnn(self,ER):
        pass

    @abstractmethod
    def FS1pp(self,ER):
        pass

    @abstractmethod
    def FS1np(self,ER):
        pass

    @abstractmethod
    def FS1pn(self,ER):
        pass

    @abstractmethod
    def FS1nn(self,ER):
        pass

    @abstractmethod
    def FS2pp(self,ER):
        pass

    @abstractmethod
    def FS2np(self,ER):
        pass

    @abstractmethod
    def FS2pn(self,ER):
        pass

    @abstractmethod
    def FS2nn(self,ER):
        pass

    @abstractmethod
    def FDpp(self,ER):
        pass

    @abstractmethod
    def FDnp(self,ER):
        pass

    @abstractmethod
    def FDpn(self,ER):
        pass

    @abstractmethod
    def FDnn(self,ER):
        pass

    @abstractmethod
    def FS1Dpp(self,ER):
        pass

    @abstractmethod
    def FS1Dnp(self,ER):
        pass

    @abstractmethod
    def FS1Dpn(self,ER):
        pass

    @abstractmethod
    def FS1Dnn(self,ER):
        pass

    ### If we wanted to we could also define the isospin basis either separately or inheriting from these guys
    ### Need to add the operators that inherit from these.... yuck. Could do it on a model by model basis, but they are nominally inherited from all of these, so fine to do here?

    def F11(self,ER,cp,cn):
        """
        Defining O1 based on form factors
        """
        return cp*cp*self.FMpp(ER)+cp*cn*self.FMpn(ER)+cn*cp*self.FMnp(ER)+cn*cn*self.FMnn(ER)
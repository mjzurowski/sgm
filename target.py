from abc import ABC, abstractmethod
import numpy as np
from constants import *

#### Fundamental properties of a target that won't change between detectors (e.g., mass, form factors)
##### Note that a lot of these functions depend on additional data. These are included in subfolders in targets/, where info.txt files can be found that give detail on their units and source

class Target(ABC):
    def __init__(self, quenching_factor=None):
        self.qf = quenching_factor

    def quenching_factor(self):
        return self.qf

    def spin_dep(self, jx):
        return (4*jx/3)*(jx+1) #### double check this defn as I think it should depend on the target
    
    @abstractmethod
    def A(self):
        """
        Mass number of target [unitless]
        """
        pass

    def Z(self):
        """
        Atomic number of target [unitless]
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

    def FDpp(self,ER):
        pass

    def FDnp(self,ER):
        pass

    def FDpn(self,ER):
        pass

    def FDnn(self,ER):
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
        rn = np.power(c1*c1+(7/3)*0.52*0.52*np.pi*np.pi-5*0.9*0.9,0.5) #fm
        q = self.Q(ER)/(197.327*1E6) # Q in units of inverse fm
        return 3*np.exp(-0.5*(0.9*q)**2)*(np.sin(q*rn)-q*rn*np.cos(q*rn))/(np.power(q*rn,3))

    ### If we wanted to we could also define the isospin basis either separately or inheriting from these guys
    ### Note that the form factors below are constructed to be called for typical NREFT DM.
    ### Some are defined with explicit g and h subscripts to match with velocity dists.
    ### Where no subscript is present should assume it is to be matched with the standard g velocity integral.

    def F11(self,ER,cp,cn):
        """
        Defining O1,1 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        p_p = cp*cp*self.FMpp(ER)
        p_n = cp*cn*self.FMpn(ER)
        n_p = cn*cp*self.FMnp(ER)
        n_n = cn*cn*self.FMnn(ER)
        return p_p+p_n+n_p+n_n
    
    def F44(self,ER,cp,cn,jx):
        """
        O4,4 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O4 operator
        """
        p_p = cp*cp*(self.FS1pp(ER)+self.FS2pp(ER))
        p_n = cp*cn*(self.FS1pn(ER)+self.FS2pn(ER))
        n_p = cn*cp*(self.FS1np(ER)+self.FS2np(ER))
        n_n = cn*cn*(self.FS1nn(ER)+self.FS2nn(ER))
        return self.spin_dep(jx)*(p_p+p_n+n_p+n_n)/16
    
    def F55(self,ER,cp,cn,jx,vm):
        """
        O5,5 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O5 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        h_p_p = cp*cp*self.FMpp(ER)
        h_p_n = cp*cn*self.FMpn(ER)
        h_n_p = cn*cp*self.FMnp(ER)
        h_n_n = cn*cn*self.FMnn(ER)
        h = self.spin_dep(jx)*np.power(self.Q(ER)/mp,2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/4

        g_p_p = cp*cp*(np.power(self.Q(ER)/mp,4.)*self.FDpp(ER)-np.power(vm*self.Q(ER)/mp,2.)*self.FMpp(ER))
        g_p_n = cp*cn*(np.power(self.Q(ER)/mp,4.)*self.FDpn(ER)-np.power(vm*self.Q(ER)/mp,2.)*self.FMpn(ER))
        g_n_p = cn*cp*(np.power(self.Q(ER)/mp,4.)*self.FDnp(ER)-np.power(vm*self.Q(ER)/mp,2.)*self.FMnp(ER))
        g_n_n = cn*cn*(np.power(self.Q(ER)/mp,4.)*self.FDnn(ER)-np.power(vm*self.Q(ER)/mp,2.)*self.FMnn(ER))
        g = self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/4

        return [g,h]
    
    def F66(self,ER,cp,cn,jx):
        """
        O6,6 operator
        cn and cp both couplings between relevant nucleon and O6 operator
        """
        p_p = cp*cp*self.FS2pp(ER)
        p_n = cp*cn*self.FS2pn(ER)
        n_p = cn*cp*self.FS2np(ER)
        n_n = cn*cn*self.FS2nn(ER)
        return self.spin_dep(jx)*np.power(self.Q(ER)/mp,4.)*(p_p+p_n+n_p+n_n)/16
    
    def F77(self,ER,cp,cn,vm):
        """
        O7,7 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O7 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        p_p = cp*cp*self.FS1pp(ER)
        p_n = cp*cn*self.FS1pn(ER)
        n_p = cn*cp*self.FS1np(ER)
        n_n = cn*cn*self.FS1nn(ER)

        h = (p_p+p_n+n_p+n_n)/8
        g = -vm*vm(p_p+p_n+n_p+n_n)/8

        return [g,h]
    
    def F10(self,ER,cp,cn,jx):
        """
        O10,10 operator
        cn and cp both couplings between relevant nucleon and O10 operator
        """
        p_p = cp*cp*self.FS2pp(ER)
        p_n = cp*cn*self.FS2pn(ER)
        n_p = cn*cp*self.FS2np(ER)
        n_n = cn*cn*self.FS2nn(ER)
        return np.power(self.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/4 ### do I need to add the converstion factor? I don't think so...
    
###### FF below here are inteference terms, and so depend on couplings to two operators
    def F45(self,ER,cp4,cn4,cp5,cn5,jx):
        """
        O4,5 operator, depends on DM spin
        """
        p_p = cp4*cp5*self.FS1Dpp(ER)
        p_n = cp4*cn5*self.FS1Dpn(ER)
        n_p = cn4*cp5*self.FS1Dnp(ER)
        n_n = cn4*cn5*self.FS1Dnn(ER)
        return 0.25*self.spin_dep(jx)*np.power(self.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)
    
    def F46(self,ER,cp4,cn4,cp6,cn6,jx):
        """
        O4,5 operator, depends on DM spin
        """
        p_p = cp4*cp6*self.FS2pp(ER)
        p_n = cp4*cn6*self.FS2pn(ER)
        n_p = cn4*cp6*self.FS2np(ER)
        n_n = cn4*cn6*self.FS2nn(ER)
        return self.spin_dep(jx)*np.power(self.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/8
    
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
        for i in range(0,len(self.eNL())):
            trans_list.append(np.interp(E_EM - self.eNL()[i], self.eTransE_E()[i][:,0], self.eTransE_E()[i][:,1]) * np.heaviside(E_EM - self.eNL()[i], 1))
        return trans_list
    
##### Photoelectric absorption
    def sigma_PE(self,E_gam):
        """
        Photoelectric absorption as a function of photon energy (E_gam) [eV]
        Output units are [cm2/atom]
        """
        pass
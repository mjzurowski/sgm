from abc import ABC, abstractmethod
import numpy as np
from constants import *
from sympy import symbols

Shell = symbols('Shell') # refers to the shell model interaction employed

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
    
    def FMpp(self,Shell,ER):
        pass

    def FMnp(self,Shell,ER):
        pass

    def FMpn(self,Shell,ER):
        pass

    def FMnn(self,Shell,ER):
        pass

    def FS1pp(self,Shell,ER):
        pass

    def FS1np(self,Shell,ER):
        pass

    def FS1pn(self,Shell,ER):
        pass

    def FS1nn(self,Shell,ER):
        pass

    def FS2pp(self,Shell,ER):
        pass

    def FS2np(self,Shell,ER):
        pass

    def FS2pn(self,Shell,ER):
        pass

    def FS2nn(self,Shell,ER):
        pass
        
    def FPhi2pp(self,Shell,ER):
        pass

    def FPhi2np(self,Shell,ER):
        pass

    def FPhi2pn(self,Shell,ER):
        pass

    def FPhi2nn(self,Shell,ER):
        pass

    def FPhipp(self,Shell,ER):
        pass

    def FPhinp(self,Shell,ER):
        pass

    def FPhipn(self,Shell,ER):
        pass

    def FPhinn(self,Shell,ER):
        pass    
    
    def FDpp(self,Shell,ER):
        pass

    def FDnp(self,Shell,ER):
        pass

    def FDpn(self,Shell,ER):
        pass

    def FDnn(self,Shell,ER):
        pass

    def FMPhi2pp(self,Shell,ER):
        pass

    def FMPhi2np(self,Shell,ER):
        pass

    def FMPhi2pn(self,Shell,ER):
        pass

    def FMPhi2nn(self,Shell,ER):
        pass
    
    def FS1Dpp(self,Shell,ER):
        pass

    def FS1Dnp(self,Shell,ER):
        pass

    def FS1Dpn(self,Shell,ER):
        pass

    def FS1Dnn(self,Shell,ER):
        pass

    def Helm(self,ER):
        c1 = 1.23*np.power(self.A(),1/3)-0.6
        s = 0.9
        rn = np.power(c1*c1+(7/3)*0.52*0.52*np.pi*np.pi-5*s*s,0.5) #fm
        q = self.Q(ER)/(197.327*1E6) # Q in units of inverse fm
        return 3*np.exp(-0.5*(s*q)**2)*(np.sin(q*rn)-q*rn*np.cos(q*rn))/(np.power(q*rn,3))

    ### If we wanted to we could also define the isospin basis either separately or inheriting from these guys
    ### Note that the form factors below are constructed to be called for typical NREFT DM.
    ### Some are defined with explicit g and h subscripts to match with velocity dists.
    ### Where no subscript is present should assume it is to be matched with the standard g velocity integral.


##############     ###############

#redefine the Fij operators based on Cirelli paper + include all of them here
# check whether cp and cn should take on different values for the cross terms, i.e. cp1, cp2, cn1, cn2
#re-write in terms of Cirelli formalism and sign convention
#check where the spin_dep(jx) terms are supposed to go
#The below are based on page 93-94 of https://arxiv.org/pdf/2104.12785
#units gives the units of the form factor component of the expressions, i.e. excluding the cp and cn coefficients. 
#The overall excpected units including the coefficients, i.e. cp cp Fij should be 

    
    def F11(self,Shell,ER,cp,cn):
        """
        Defining O1,1 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        p_p = cp*cp*self.FMpp(Shell,ER)
        p_n = cp*cn*self.FMpn(Shell,ER)
        n_p = cn*cp*self.FMnp(Shell,ER)
        n_n = cn*cn*self.FMnn(Shell,ER)
        return p_p+p_n+n_p+n_n   #units=untiless
        
    def F33(self,Shell,ER,cp,cn):
        """
        Defining O3,3 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        
        h_p_p = cp*cp*self.FS1pp(Shell,ER)
        h_p_n = cp*cn*self.FS1pn(Shell,ER)
        h_n_p = cn*cp*self.FS1np(Shell,ER)
        h_n_n = cn*cn*self.FS1nn(Shell,ER)
        h = self.spin_dep(jx)*np.power(self.Q(ER),2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/8  #units=ev^2

        g_p_p = cp*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2pp(Shell,ER)-np.power(vm*self.Q(ER),2.)*self.FS1pp(Shell,ER)/2)
        g_p_n = cp*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2pn(Shell,ER)-np.power(vm*self.Q(ER),2.)*self.FS1pn(Shell,ER)/2)
        g_n_p = cn*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2np(Shell,ER)-np.power(vm*self.Q(ER),2.)*self.FS1np(Shell,ER)/2)
        g_n_n = cn*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2nn(Shell,ER)-np.power(vm*self.Q(ER),2.)*self.FS1nn(Shell,ER)/2)
        g = self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/4  #units=ev^2

        return [g,h]


    
    def F44(self,Shell,ER,cp,cn,jx):
        """
        O4,4 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O4 operator
        """
        p_p = cp*cp*(self.FS1pp(ER)+self.FS2pp(ER))
        p_n = cp*cn*(self.FS1pn(ER)+self.FS2pn(ER))
        n_p = cn*cp*(self.FS1np(ER)+self.FS2np(ER))
        n_n = cn*cn*(self.FS1nn(ER)+self.FS2nn(ER))
        return self.spin_dep(jx)*(p_p+p_n+n_p+n_n)/16  #units= unitless

        
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
        h = self.spin_dep(jx)*np.power(self.Q(ER),2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/4  #units=ev^2

        g_p_p = cp*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.FDpp(ER)-np.power(vm*self.Q(ER),2.)*self.FMpp(ER))
        g_p_n = cp*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.FDpn(ER)-np.power(vm*self.Q(ER),2.)*self.FMpn(ER))
        g_n_p = cn*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.FDnp(ER)-np.power(vm*self.Q(ER),2.)*self.FMnp(ER))
        g_n_n = cn*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.FDnn(ER)-np.power(vm*self.Q(ER),2.)*self.FMnn(ER))
        g = self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/4   #units= ev^2

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
        return self.spin_dep(jx)*np.power(self.Q(ER),4.)*(p_p+p_n+n_p+n_n)/16  #units=ev^4

    
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

        h = (p_p+p_n+n_p+n_n)/8   #units=unitless
        g = -vm*vm(p_p+p_n+n_p+n_n)/8  #units=unitless

        return [g,h]

    
    def F88(self,ER,cp,cn,jx,vm):
        """
        O8,8 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O5 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        #there was a q^2 included in h that I don't think is supposed to be there, I've removed it
        
        h_p_p = cp*cp*self.FMpp(ER)
        h_p_n = cp*cn*self.FMpn(ER)
        h_n_p = cn*cp*self.FMnp(ER)
        h_n_n = cn*cn*self.FMnn(ER)
        h = self.spin_dep(jx)*(h_p_p+h_p_n+h_n_p+h_n_n)/4  #units=unitless

        g_p_p = cp*cp*(np.power(self.Q(ER)/mp,2.)*self.FDpp(ER)-np.power(vm,2.)*self.FMpp(ER))
        g_p_n = cp*cn*(np.power(self.Q(ER)/mp,2.)*self.FDpn(ER)-np.power(vm,2.)*self.FMpn(ER))
        g_n_p = cn*cp*(np.power(self.Q(ER)/mp,2.)*self.FDnp(ER)-np.power(vm,2.)*self.FMnp(ER))
        g_n_n = cn*cn*(np.power(self.Q(ER)/mp,2.)*self.FDnn(ER)-np.power(vm,2.)*self.FMnn(ER))
        g = self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/4  #units=unitless

        return [g,h]

    def F99(self,ER,cp,cn,jx):
        """
        O9,9 operator
        cn and cp both couplings between relevant nucleon and O6 operator
        """
        p_p = cp*cp*self.FS1pp(ER)
        p_n = cp*cn*self.FS1pn(ER)
        n_p = cn*cp*self.FS1np(ER)
        n_n = cn*cn*self.FS1nn(ER)
        return self.spin_dep(jx)*np.power(self.Q(ER),2.)*(p_p+p_n+n_p+n_n)/16  #units=ev^2
    
    
    def F1010(self,ER,cp,cn,jx):
        """
        O10,10 operator
        cn and cp both couplings between relevant nucleon and O10 operator
        """
        p_p = cp*cp*self.FS2pp(ER)
        p_n = cp*cn*self.FS2pn(ER)
        n_p = cn*cp*self.FS2np(ER)
        n_n = cn*cn*self.FS2nn(ER)
        return np.power(self.Q(ER),2.)*(p_p+p_n+n_p+n_n)/4  #units=ev^2
        
    ### do I need to add the converstion factor? I don't think so...

    def F1111(self,ER,cp,cn):
        """
        Defining O11,11 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        p_p = cp*cp*self.FMpp(ER)
        p_n = cp*cn*self.FMpn(ER)
        n_p = cn*cp*self.FMnp(ER)
        n_n = cn*cn*self.FMnn(ER)
        return np.power(self.Q(ER),2.)*(p_p+p_n+n_p+n_n)/4  #units=ev^2

    
    def F1212(self,ER,cp,cn):
        """
        Defining O12,12 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        
        h_p_p = cp*cp*(self.FS1pp(ER)/2 + self.FS2pp(ER))
        h_p_n = cp*cn*(self.FS1pn(ER)/2 + self.FS2pn(ER))
        h_n_p = cn*cp*(self.FS1np(ER)/2 + self.FS2np(ER))
        h_n_n = cn*cn*(self.FS1nn(ER)/2 + self.FS2nn(ER))
        h = self.spin_dep(jx)*(h_p_p+h_p_n+h_n_p+h_n_n)/16  #units=unitless

        g_p_p = cp*cp*(np.power(self.Q(ER)/mp,2.)*(self.Phi2pp(ER) + self.Phipp(ER))-np.power(vm,2.)*(self.FS1pp(ER)/2 + self.FS2pp(ER)))
        g_p_n = cp*cn*(np.power(self.Q(ER)/mp,2.)*(self.Phi2pn(ER) + self.Phipn(ER))-np.power(vm,2.)*(self.FS1pn(ER)/2 + self.FS2pn(ER)))
        g_n_p = cn*cp*(np.power(self.Q(ER)/mp,2.)*(self.Phi2np(ER) + self.Phinp(ER))-np.power(vm,2.)*(self.FS1np(ER)/2 + self.FS2np(ER)))
        g_n_n = cn*cn*(np.power(self.Q(ER)/mp,2.)*(self.Phi2nn(ER) + self.Phinn(ER))-np.power(vm,2.)*(self.FS1nn(ER)/2 + self.FS2nn(ER)))
        g = self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  #units=untiless

        return [g,h]
   
    
    def F1313(self,ER,cp,cn):
        """
        Defining O13,13 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        
        h_p_p = cp*cp*self.FS2pp(ER)
        h_p_n = cp*cn*self.FS2pn(ER)
        h_n_p = cn*cp*self.FS2np(ER)
        h_n_n = cn*cn*self.FS2nn(ER)
        h = self.spin_dep(jx)*np.power(self.Q(ER),2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/16  #units=ev^2

        g_p_p = cp*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phipp(ER)-np.power(vm*self.Q(ER),2.)*self.FS2pp(ER))
        g_p_n = cp*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phipn(ER)-np.power(vm*self.Q(ER),2.)*self.FS2pn(ER))
        g_n_p = cn*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phinp(ER)-np.power(vm*self.Q(ER),2.)*self.FS2np(ER))
        g_n_n = cn*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phinn(ER)-np.power(vm*self.Q(ER),2.)*self.FS2nn(ER))
        g = self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  #units=ev^2

        return [g,h]


    def F1414(self,ER,cp,cn,vm):
        """
        O14,14 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O7 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        p_p = cp*cp*self.FS1pp(ER)
        p_n = cp*cn*self.FS1pn(ER)
        n_p = cn*cp*self.FS1np(ER)
        n_n = cn*cn*self.FS1nn(ER)

        h = np.power(self.Q(ER),2.)*(p_p+p_n+n_p+n_n)/32    #units=ev^2
        g = -vm*vm*np.power(self.Q(ER),2.)*(p_p+p_n+n_p+n_n)/32   #units=ev^2

        return [g,h]
    

    def F1515(self,ER,cp,cn):
        """
        Defining O15,15 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        
        h_p_p = cp*cp*self.FS1pp(ER)
        h_p_n = cp*cn*self.FS1pn(ER)
        h_n_p = cn*cp*self.FS1np(ER)
        h_n_n = cn*cn*self.FS1nn(ER)
        h = self.spin_dep(jx)*np.power(self.Q(ER),4.)*(h_p_p+h_p_n+h_n_p+h_n_n)/32  #units=ev^4

        g_p_p = cp*cp*(np.power(self.Q(ER)/(mp**(1/3)),6.)*self.Phi2pp(ER)-np.power(math.sqrt(vm)*self.Q(ER),4.)*self.FS1pp(ER)/2)
        g_p_n = cp*cn*(np.power(self.Q(ER)/(mp**(1/3)),6.)*self.Phi2pn(ER)-np.power(math.sqrt(vm)*self.Q(ER),4.)*self.FS1pn(ER)/2)
        g_n_p = cn*cp*(np.power(self.Q(ER)/(mp**(1/3)),6.)*self.Phi2np(ER)-np.power(math.sqrt(vm)*self.Q(ER),4.)*self.FS1np(ER)/2)
        g_n_n = cn*cn*(np.power(self.Q(ER)/(mp**(1/3)),6.)*self.Phi2nn(ER)-np.power(math.sqrt(vm)*self.Q(ER),4.)*self.FS1nn(ER)/2)
        g = self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  #units=ev^4

        return [g,h]    

    
    
###### FF below here are inteference terms, and so depend on couplings to two operators
    def F13(self,ER,cp1,cn1,cp3,cn3,jx):
        """
        O1,3 operator, depends on DM spin
        """
        p_p = cp1*cp3*self.FMPhi2pp(ER)
        p_n = cp1*cn3*self.FMPhi2pn(ER)
        n_p = cn1*cp3*self.FMPhi2np(ER)
        n_n = cn1*cn3*self.FMPhi2nn(ER)
        return -self.spin_dep(jx)*np.power(self.Q(ER)/math.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/2   #units=ev
    

    def F45(self,ER,cp4,cn4,cp5,cn5,jx):
        """
        O4,5 operator, depends on DM spin
        """
        p_p = cp4*cp5*self.FS1Dpp(ER)
        p_n = cp4*cn5*self.FS1Dpn(ER)
        n_p = cn4*cp5*self.FS1Dnp(ER)
        n_n = cn4*cn5*self.FS1Dnn(ER)
        return -self.spin_dep(jx)*np.power(self.Q(ER)/math.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8   #units=ev
    
    def F46(self,ER,cp4,cn4,cp6,cn6,jx):
        """
        O4,6 operator, depends on DM spin
        """
        p_p = cp4*cp6*self.FS2pp(ER)
        p_n = cp4*cn6*self.FS2pn(ER)
        n_p = cn4*cp6*self.FS2np(ER)
        n_n = cn4*cn6*self.FS2nn(ER)
        return self.spin_dep(jx)*np.power(self.Q(ER),2.)*(p_p+p_n+n_p+n_n)/16  #units=ev^2

#double check 9,8 and the order of the S1 D form factor
    def F98(self,ER,cp9,cn9,cp8,cn8,jx):
        """
        O9,8 operator, depends on DM spin
        """
        p_p = cp9*cp8*self.FS1Dpp(ER)
        p_n = cp9*cn8*self.FS1Dpn(ER)
        n_p = cn9*cp8*self.FS1Dnp(ER)
        n_n = cn9*cn8*self.FS1Dnn(ER)
        return self.spin_dep(jx)*np.power(self.Q(ER)/math.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8   #units=ev


    def F1112(self,ER,cp11,cn11,cp12,cn12,jx):
        """
        O11,12 operator, depends on DM spin
        """
        p_p = cp11*cp12*self.FMPhi2pp(ER)
        p_n = cp11*cn12*self.FMPhi2pn(ER)
        n_p = cn11*cp12*self.FMPhi2np(ER)
        n_n = cn11*cn12*self.FMPhi2nn(ER)
        return -self.spin_dep(jx)*np.power(self.Q(ER)/math.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8   #units=ev


    def F1115(self,ER,cp11,cn11,cp15,cn15,jx):
        """
        O11,15 operator, depends on DM spin
        """
        p_p = cp11*cp15*self.FMPhi2pp(ER)
        p_n = cp11*cn15*self.FMPhi2pn(ER)
        n_p = cn11*cp15*self.FMPhi2np(ER)
        n_n = cn11*cn15*self.FMPhi2nn(ER)
        return self.spin_dep(jx)*np.power(self.Q(ER)/(mp**(1/4)),4.)*(p_p+p_n+n_p+n_n)/8   #units=ev^3


    def F1215(self,ER,cp12,cn12,cp15,cn15,jx):
        """
        Defining O12,15 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        
        h_p_p = cp*cp*self.FS1pp(ER)
        h_p_n = cp*cn*self.FS1pn(ER)
        h_n_p = cn*cp*self.FS1np(ER)
        h_n_n = cn*cn*self.FS1nn(ER)
        h = -self.spin_dep(jx)*np.power(self.Q(ER),2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/32  #units=ev^2

        g_p_p = cp*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2pp(ER)-np.power(vm*self.Q(ER),2.)*self.FS1pp(ER)/2)
        g_p_n = cp*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2pn(ER)-np.power(vm*self.Q(ER),2.)*self.FS1pn(ER)/2)
        g_n_p = cn*cp*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2np(ER)-np.power(vm*self.Q(ER),2.)*self.FS1np(ER)/2)
        g_n_n = cn*cn*(np.power(self.Q(ER)/math.sqrt(mp),4.)*self.Phi2nn(ER)-np.power(vm*self.Q(ER),2.)*self.FS1nn(ER)/2)
        g = -self.spin_dep(jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  #units=ev^2

        return [g,h]    


#####################


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




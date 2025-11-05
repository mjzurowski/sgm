from abc import ABC, abstractmethod
import numpy as np
from constants import *
from models.couplings import *

"""
Class defintions for all the DM form factors, following the Anand formalism 
(i.e., unitless form factors as per arxiv 1308.6288)
These models are for use where you want to plot/constrain an experimental cross section
"""

"""
Create an abstract class to automatically calc the appropriate cross sections etc for an elastic DM model
This is formatted so you can also use it in combination with other inelastic models (e.g., inelastic DM, Migdal)
"""
class Anand(ABC):
    def vmin(self,Target,mX,ER,**kwargs):
        """
        Minimum velocity that can produce recoil energy ER in [km/s]
        [Target]: target nucleus
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: [km/s]
        """
        return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)

    @abstractmethod
    def FF(self, Target, ER, vm):
        """
        Target = Target type object
        [ER] = [eV] DM recoil energy
        [vm] = minimum velocity [km/s] (if needed)

        Output units: unitless
        Note that if this operator depends on velocity this should return a list of the g and h form factors
        """
        pass
    
    def dsigdER(self,Target,ER,mX,sig,vm):
        """
        Differential cross section in units of [cm^2]/[eV]
        """
        FF = self.FF(Target,ER,vm/kms)
        cross_sec = sig 
        try:
            _dsigdER0 = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            _dsigdER1 = cross_sec*FF[1]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return [_dsigdER0,_dsigdER1]
        except:
            _dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return _dsigdER
            

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
        vm = self.vmin(Target,mX,ER)
        _dsigdER = self.dsigdER(Target,ER,mX,sig,vm)
        try: 
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(_dsigdER[0]*VelDist.gdist(vm) + _dsigdER[1]*VelDist.hdist(vm)) 
        except:
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*VelDist.gdist(vm)*_dsigdER
                

########################################################
class AnandF1(Anand):
    def __init__(self, cp, cn,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn

    def FF(self, Target, ER, vm):
        """
        Form factor expression for O1
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp*self.cp*Target.FMpp(ER)
            p_n = self.cp*self.cn*Target.FMpn(ER)
            n_p = self.cn*self.cp*Target.FMnp(ER)
            n_n = self.cn*self.cn*Target.FMnn(ER)
            return p_p+p_n+n_p+n_n


class AnandF3(Anand):
    def __init__(self, cp, cn,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn

    def FF(self, Target, ER,vm):
        """
        Form factor expression for O3
        Target = Target type object
        [ER] = [eV] DM recoil energy
        [vm] = [unitless] minimun velocity for recoil of ER

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            h_p_p = self.cp*self.cp*Target.FS1pp(ER)
            h_p_n = self.cp*self.cn*Target.FS1pn(ER)
            h_n_p = self.cn*self.cp*Target.FS1np(ER)
            h_n_n = self.cn*self.cn*Target.FS1nn(ER)
            h = np.power(Target.Q(ER)/mp,2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/8  #unitless

            g_p_p = self.cp*self.cp*(np.power(Target.Q(ER)/mp,4.)*Target.FPhi2pp(ER)/4-np.power(vm*Target.Q(ER)/mp,2.)*Target.FS1pp(ER)/8)
            g_p_n = self.cp*self.cn*(np.power(Target.Q(ER)/mp,4.)*Target.FPhi2pn(ER)/4-np.power(vm*Target.Q(ER)/mp,2.)*Target.FS1pn(ER)/8)
            g_n_p = self.cn*self.cp*(np.power(Target.Q(ER)/mp,4.)*Target.FPhi2np(ER)/4-np.power(vm*Target.Q(ER)/mp,2.)*Target.FS1np(ER)/8)
            g_n_n = self.cn*self.cn*(np.power(Target.Q(ER)/mp,4.)*Target.FPhi2nn(ER)/4-np.power(vm*Target.Q(ER)/mp,2.)*Target.FS1nn(ER)/8)
            g = (g_p_p+g_p_n+g_n_p+g_n_n)  #unitless

            return [g,h]
        
class AnandF4(Anand):
    def __init__(self, cp, cn,jx,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn
        self.jx = jx


    def FF(self, Target, ER, vm):
        """
        Form factor expression for O4
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp*self.cp*(Target.FS1pp(ER)+Target.FS2pp(ER))
            p_n = self.cp*self.cn*(Target.FS1pn(ER)+Target.FS2pn(ER))
            n_p = self.cn*self.cp*(Target.FS1np(ER)+Target.FS2np(ER))
            n_n = self.cn*self.cn*(Target.FS1nn(ER)+Target.FS2nn(ER))

            return Target.spin_dep(self.jx)*(p_p+p_n+n_p+n_n)/16

        
class AnandF5(Anand):
    def __init__(self, cp, cn, jx,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn
        self.jx = jx

    def FF(self, Target, ER, vm):
        """
        Form factor expression for O5
        Target = Target type object
        [ER] = [eV] DM recoil energy
        [vm] = [unitless] minimun velocity for recoil of ER

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            h_p_p = self.cp*self.cp*Target.FMpp(ER)
            h_p_n = self.cp*self.cn*Target.FMpn(ER)
            h_n_p = self.cn*self.cp*Target.FMnp(ER)
            h_n_n = self.cn*self.cn*Target.FMnn(ER)
            h = Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/4

            g_p_p = self.cp*self.cp*(np.power(Target.Q(ER)/mp,4.)*Target.FDpp(ER)-np.power(vm*Target.Q(ER)/mp,2.)*Target.FMpp(ER))
            g_p_n = self.cp*self.cn*(np.power(Target.Q(ER)/mp,4.)*Target.FDpn(ER)-np.power(vm*Target.Q(ER)/mp,2.)*Target.FMpn(ER))
            g_n_p = self.cn*self.cp*(np.power(Target.Q(ER)/mp,4.)*Target.FDnp(ER)-np.power(vm*Target.Q(ER)/mp,2.)*Target.FMnp(ER))
            g_n_n = self.cn*self.cn*(np.power(Target.Q(ER)/mp,4.)*Target.FDnn(ER)-np.power(vm*Target.Q(ER)/mp,2.)*Target.FMnn(ER))
            g = Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/4

            return [g,h]
        
class AnandF6(Anand):
    def __init__(self, cp, cn, jx,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn
        self.jx = jx

    def FF(self, Target, ER, vm):
        """
        Form factor expression for O6
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp*self.cp*(Target.FS2pp(ER))
            p_n = self.cp*self.cn*(Target.FS2pn(ER))
            n_p = self.cn*self.cp*(Target.FS2np(ER))
            n_n = self.cn*self.cn*(Target.FS2nn(ER))

            return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,4.)*(p_p+p_n+n_p+n_n)/16
        
class AnandF7(Anand):
    def __init__(self, cp, cn,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn

    def FF(self,Target,ER,vm):
        """
        Form factor expression for O7
        Target = Target type object
        [ER] = [eV] DM recoil energy
        [vm] = [unitless] minimun velocity for recoil of ER

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp*self.cp*Target.FS1pp(ER)
            p_n = self.cp*self.cn*Target.FS1pn(ER)
            n_p = self.cn*self.cp*Target.FS1np(ER)
            n_n = self.cn*self.cn*Target.FS1nn(ER)

            h = (p_p+p_n+n_p+n_n)/8
            g = -vm*vm*(p_p+p_n+n_p+n_n)/8

            return [g,h]
        
class AnandF8(Anand):
    def __init__(self, cp, cn, jx,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn
        self.jx = jx

    def FF(self,Target,ER,vm):
        """
        Form factor expression for O8
        Target = Target type object
        [ER] = [eV] DM recoil energy
        [vm] = [unitless] minimun velocity for recoil of ER

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            g_p_p = self.cp*self.cp*(pow(Target.Q(ER)/mp,2)*Target.FDpp(ER)-vm*vm*Target.FMpp(ER))
            g_p_n = self.cp*self.cn*(pow(Target.Q(ER)/mp,2)*Target.FDpn(ER)-vm*vm*Target.FMpn(ER))
            g_n_p = self.cn*self.cp*(pow(Target.Q(ER)/mp,2)*Target.FDnp(ER)-vm*vm*Target.FMnp(ER))
            g_n_n = self.cn*self.cn*(pow(Target.Q(ER)/mp,2)*Target.FDnn(ER)-vm*vm*Target.FMnn(ER))
            g = Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/4
            
            h_p_p = self.cp*self.cp*Target.FMpp(ER)
            h_p_n = self.cp*self.cn*Target.FMpn(ER)
            h_n_p = self.cn*self.cp*Target.FMnp(ER)
            h_n_n = self.cn*self.cn*Target.FMnn(ER)
            h = Target.spin_dep(self.jx)*(h_p_p+h_p_n+h_n_p+h_n_n)/4

            return [g,h]
    
        
class AnandF9(Anand):
    def __init__(self, cp, cn, jx,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn
        self.jx = jx

    def FF(self,Target,ER,vm):
        """
        Form factor expression for O9
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp*self.cp*Target.FS1pp(ER)
            p_n = self.cp*self.cn*Target.FS1pn(ER)
            n_p = self.cn*self.cp*Target.FS1np(ER)
            n_n = self.cn*self.cn*Target.FS1nn(ER)

            return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/16 

        
class AnandF10(Anand):
    def __init__(self, cp, cn,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn

    def FF(self,Target,ER,vm):
        """
        Form factor expression for O10
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp*self.cp*Target.FS2pp(ER)
            p_n = self.cp*self.cn*Target.FS2pn(ER)
            n_p = self.cn*self.cp*Target.FS2np(ER)
            n_n = self.cn*self.cn*Target.FS2nn(ER)

            return np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/4 

class AnandF11(Anand):
    def __init__(self, cp, cn,jx,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn
        self.jx = jx

    def FF(self, Target, ER,vm):
        """
        Form factor expression for O11
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp*self.cp*Target.FMpp(ER)
            p_n = self.cp*self.cn*Target.FMpn(ER)
            n_p = self.cn*self.cp*Target.FMnp(ER)
            n_n = self.cn*self.cn*Target.FMnn(ER)
            return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/4
        
class AnandF12(Anand):
    def __init__(self, cp, cn, jx,norm="p"):
        """
        Initialise with a set of cp and cn values
        cp: coupling to proton
        cn: coupling to neutron
        norm: normalisation method used (default is to let cp==1). This will impact the physical meaning of your cross section.
            (if you aren't sure what that means, just use "p", which is the default for usual direct detection sensitivity plots)
        """
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp = cp/cp
            self.cn = cn/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp = cp/cn
            self.cn = cn/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp = cp/np.sqrt(cp*cp+cn*cn)
            self.cn = cn/np.sqrt(cp*cp+cn*cn)
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp = cp
            self.cn = cn
        self.jx = jx

    def FF(self,Target,ER,vm): 
        """
        Form factor expression for O12
        Target = Target type object
        [ER] = [eV] DM recoil energy
        [vm] = [unitless] minimun velocity for recoil of ER

        Output units: unitless
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            h_p_p = self.cp*self.cp*(Target.FS1pp(ER)/2 + Target.FS2pp(ER))
            h_p_n = self.cp*self.cn*(Target.FS1pn(ER)/2 + Target.FS2pn(ER))
            h_n_p = self.cn*self.cp*(Target.FS1np(ER)/2 + Target.FS2np(ER))
            h_n_n = self.cn*self.cn*(Target.FS1nn(ER)/2 + Target.FS2nn(ER))
            h = Target.spin_dep(self.jx)*(h_p_p+h_p_n+h_n_p+h_n_n)/16 

            g_p_p = self.cp*self.cp*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2pp(ER) + Target.FPhipp(ER))-np.power(vm,2.)*(Target.FS1pp(ER)/2 + Target.FS2pp(ER)))
            g_p_n = self.cp*self.cn*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2pn(ER) + Target.FPhipn(ER))-np.power(vm,2.)*(Target.FS1pn(ER)/2 + Target.FS2pn(ER)))
            g_n_p = self.cn*self.cp*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2np(ER) + Target.FPhinp(ER))-np.power(vm,2.)*(Target.FS1np(ER)/2 + Target.FS2np(ER)))
            g_n_n = self.cn*self.cn*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2nn(ER) + Target.FPhinn(ER))-np.power(vm,2.)*(Target.FS1nn(ER)/2 + Target.FS2nn(ER)))
            g = Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  

            return [g,h]
        
###### Need to double check units, signs, and normalisation for these       
class AnandF1F3(Anand):
    def __init__(self, cp1, cn1, cp3, cn3, norm = "p"):
        c0 = np.sqrt(cp1**2 + cn1**2 + cp3**2 + cn3**2)
        cp = np.sqrt(cp1**2 + cp3**2)
        cn = np.sqrt(cn1**2 + cn3**2)
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp1 = cp1/cp
            self.cn1 = cn1/cp
            self.cp3 = cp3/cp
            self.cn3 = cn3/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp1 = cp1/cn
            self.cn1 = cn1/cn
            self.cp3 = cp3/cn
            self.cn3 = cn3/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp1 = cp1/c0
            self.cn1 = cn1/c0
            self.cp3 = cp3/c0
            self.cn3 = cn3/c0
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp1 = cp1
            self.cn1 = cn1
            self.cp3 = cp3
            self.cn3 = cn3
        self.jx = jx

    def FF(self, Target, ER,vm):
        """
        Form factor expression for interference of O1 and O3
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn1==self.cp1==self.cn3==self.cp3==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp1*self.cp3*Target.FMPhi2pp(ER)
            p_n = self.cp1*self.cn3*Target.FMPhi2pn(ER)
            n_p = self.cn1*self.cp3*Target.FMPhi2np(ER)
            n_n = self.cn1*self.cn3*Target.FMPhi2nn(ER)
            return 0.5*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)   


class AnandF4F5(Anand):
    def __init__(self, cp4, cn4, cp5, cn5, jx, norm = "p"):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        c0 = np.sqrt(cp4**2 + cn4**2 + cp5**2 + cn5**2)
        cp = np.sqrt(cp4**2 + cp5**2)
        cn = np.sqrt(cn4**2 + cn5**2)
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp4 = cp4/cp
            self.cn4 = cn4/cp
            self.cp5 = cp5/cp
            self.cn5 = cn5/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp4 = cp4/cn
            self.cn4 = cn4/cn
            self.cp5 = cp5/cn
            self.cn5 = cn5/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp4 = cp4/c0
            self.cn4 = cn4/c0
            self.cp5 = cp5/c0
            self.cn5 = cn5/c0
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp4 = cp4
            self.cn4 = cn4
            self.cp5 = cp5
            self.cn5 = cn5
        self.jx = jx

    def FF(self, Target, ER,vm):
        """
        Form factor expression for interference of O4 and O5
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn4==self.cp4==self.cn5==self.cp5==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp4*self.cp5*Target.FS1Dpp(ER)
            p_n = self.cp4*self.cn5*Target.FS1Dpn(ER)
            n_p = self.cn4*self.cp5*Target.FS1Dnp(ER)
            n_n = self.cn4*self.cn5*Target.FS1Dnn(ER)
            return 0.5*0.25*Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)
        

class AnandF4F6(Anand):
    def __init__(self, cp4, cn4, cp6, cn6, jx, norm = "p"):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        c0 = np.sqrt(cp4**2 + cn4**2 + cp6**2 + cn6**2)
        cp = np.sqrt(cp4**2 + cp6**2)
        cn = np.sqrt(cn4**2 + cn6**2)
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp4 = cp4/cp
            self.cn4 = cn4/cp
            self.cp6 = cp6/cp
            self.cn6 = cn6/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp4 = cp4/cn
            self.cn4 = cn4/cn
            self.cp6 = cp6/cn
            self.cn6 = cn6/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp4 = cp4/c0
            self.cn4 = cn4/c0
            self.cp6 = cp6/c0
            self.cn6 = cn6/c0
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp4 = cp4
            self.cn4 = cn4
            self.cp6 = cp6
            self.cn6 = cn6
        self.jx = jx

    def FF(self, Target, ER,vm):
        """
        Form factor expression for interference of O4 and O5
        Target = Target type object
        [ER] = [eV] DM recoil energy

        Output units: unitless
        """
        if(self.cn4==self.cp4==self.cn6==self.cp6==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp4*self.cp6*Target.FS2pp(ER)
            p_n = self.cp4*self.cn6*Target.FS2pn(ER)
            n_p = self.cn4*self.cp6*Target.FS2np(ER)
            n_n = self.cn4*self.cn6*Target.FS2nn(ER)
            return 0.5*Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/8

class AnandF8F9(Anand):
    def __init__(self, cp8, cn8, cp9, cn9, jx, norm = "p"):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        
        c0 = np.sqrt(cp8**2 + cn8**2 + cp9**2 + cn9**2)
        cp = np.sqrt(cp8**2 + cp9**2)
        cn = np.sqrt(cn8**2 + cn9**2)
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp8 = cp8/cp
            self.cn8 = cn8/cp
            self.cp9 = cp9/cp
            self.cn9 = cn9/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp8 = cp8/cn
            self.cn8 = cn8/cn
            self.cp9 = cp9/cn
            self.cn9 = cn9/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp8 = cp8/c0
            self.cn8 = cn8/c0
            self.cp9 = cp9/c0
            self.cn9 = cn9/c0
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp8 = cp8
            self.cn8 = cn8
            self.cp9 = cp9
            self.cn9 = cn9
        self.jx = jx

    
    def FF(self, Target, ER,vm):
        """
        Form factor expression for interference of O8 and O9
        """
        if(self.cn8==self.cp8==self.cn9==self.cp9==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp8*self.cp9*Target.FS1Dpp(ER)
            p_n = self.cp8*self.cn9*Target.FS1Dnp(ER)
            n_p = self.cn8*self.cp9*Target.FS1Dpn(ER)
            n_n = self.cn8*self.cn9*Target.FS1Dnn(ER)
            return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV

class AnandF9F8(Anand):
    def __init__(self, cp9, cn9, cp8, cn8, jx, norm = "p"):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        
        c0 = np.sqrt(cp8**2 + cn8**2 + cp9**2 + cn9**2)
        cp = np.sqrt(cp8**2 + cp9**2)
        cn = np.sqrt(cn8**2 + cn9**2)
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp8 = cp8/cp
            self.cn8 = cn8/cp
            self.cp9 = cp9/cp
            self.cn9 = cn9/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp8 = cp8/cn
            self.cn8 = cn8/cn
            self.cp9 = cp9/cn
            self.cn9 = cn9/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp8 = cp8/c0
            self.cn8 = cn8/c0
            self.cp9 = cp9/c0
            self.cn9 = cn9/c0
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp8 = cp8
            self.cn8 = cn8
            self.cp9 = cp9
            self.cn9 = cn9
        self.jx = jx

    
    def FF(self, Target, ER,vm):
        """
        Form factor expression for interference of O8 and O9
        """
        if(self.cn9==self.cp9==self.cn8==self.cp8==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp9*self.cp8*Target.FS1Dpp(ER)
            p_n = self.cp9*self.cn8*Target.FS1Dpn(ER)
            n_p = self.cn9*self.cp8*Target.FS1Dnp(ER)
            n_n = self.cn9*self.cn8*Target.FS1Dnn(ER)        
            return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV


class AnandF11F12(Anand):
    def __init__(self, cp11, cn11, cp12, cn12, jx, norm = "p"):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        c0 = np.sqrt(cp11**2 + cn11**2 + cp12**2 + cn12**2)
        cp = np.sqrt(cp11**2 + cp12**2)
        cn = np.sqrt(cn11**2 + cn12**2)
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp11 = cp11/cp
            self.cn11 = cn11/cp
            self.cp12 = cp12/cp
            self.cn12 = cn12/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp11 = cp11/cn
            self.cn11 = cn11/cn
            self.cp12 = cp12/cn
            self.cn12 = cn12/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp11 = cp11/c0
            self.cn11 = cn11/c0
            self.cp12 = cp12/c0
            self.cn12 = cn12/c0
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp11 = cp11
            self.cn11 = cn11
            self.cp12 = cp12
            self.cn12 = cn12
        self.jx = jx

    
    def FF(self, Target, ER,vm):
        """
        Form factor expression for interference of O11 and O12
        """
        if(self.cn11==self.cp11==self.cn12==self.cp12==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp11*self.cp12*Target.FMPhi2pp(ER)
            p_n = self.cp11*self.cn12*Target.FMPhi2pn(ER)
            n_p = self.cn11*self.cp12*Target.FMPhi2np(ER)
            n_n = self.cn11*self.cn12*Target.FMPhi2nn(ER)
            return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV


class AnandF12F11(Anand):
    def __init__(self, cp12, cn12, cp11, cn11, jx, norm = "p"):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        c0 = np.sqrt(cp11**2 + cn11**2 + cp12**2 + cn12**2)
        cp = np.sqrt(cp11**2 + cp12**2)
        cn = np.sqrt(cn11**2 + cn12**2)
        if norm=="p":
            print("Normalising wrt the proton coupling. Make sure to use proton cross section for EFT matching")
            self.cp11 = cp11/cp
            self.cn11 = cn11/cp
            self.cp12 = cp12/cp
            self.cn12 = cn12/cp
        if norm=="n":
            print("Normalising wrt the neutron coupling. Make sure to use neutron cross section for EFT matching")
            self.cp11 = cp11/cn
            self.cn11 = cn11/cn
            self.cp12 = cp12/cn
            self.cn12 = cn12/cn
        if norm=="vec":
            print("Normalising wrt the nucleon 'vector'. Make sure to use nucleon vector cross section for EFT matching")
            self.cp11 = cp11/c0
            self.cn11 = cn11/c0
            self.cp12 = cp12/c0
            self.cn12 = cn12/c0
        if norm=="none":
            print("No normalisation will be applied to the couplings. I hope you've carefully normalised them yourself...")
            self.cp11 = cp11
            self.cn11 = cn11
            self.cp12 = cp12
            self.cn12 = cn12
        self.jx = jx

    
    def FF(self, Target, ER,vm):
        """
        Form factor expression for interference of O11 and O12
        """
        if(self.cn11==self.cp11==self.cn12==self.cp12==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            p_p = self.cp12*self.cp11*Target.FMPhi2pp(ER)
            p_n = self.cp12*self.cn11*Target.FMPhi2np(ER)
            n_p = self.cn12*self.cp11*Target.FMPhi2pn(ER)
            n_n = self.cn12*self.cn11*Target.FMPhi2nn(ER)
            return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV
      

#### Need to add function to call the appropriate FFs based on some high energy coupling.
class AnandFull:
    """
    This is a class in and of itself to account for the structure required to combine all the FFs
    """
    def __init__(self, cq, jx,mX,Lam,norm="p"):
        # get couplings
        c1 = c1_NR(cq, mX, Lam, norm)
        c4 = c4_NR(cq, mX, Lam, norm)
        c6 = c6_NR(cq, mX, Lam, norm)
        c7 = c7_NR(cq, mX, Lam, norm)
        c8 = c8_NR(cq, mX, Lam, norm)
        c9 = c9_NR(cq, mX, Lam, norm)
        c10 = c10_NR(cq, mX, Lam, norm)
        c11 = c11_NR(cq, mX, Lam, norm)
        c12 = c12_NR(cq, mX, Lam, norm)
        self.mX = mX
        self.sig = sigma_from_EFT(cq,mX,Lam,norm)

        # get form factor model objects
        self.F1 = AnandF1(c1[0],c1[1],norm="none")
        self.F4 = AnandF4(c4[0],c4[1],jx,norm="none")
        self.F6 = AnandF6(c6[0],c6[1],jx,norm="none")
        self.F7 = AnandF7(c7[0],c7[1],norm="none")
        self.F8 = AnandF8(c8[0],c8[1],jx,norm="none")
        self.F9 = AnandF9(c9[0],c9[1],jx,norm="none")
        self.F10 = AnandF10(c10[0],c10[1],norm="none")
        self.F11 = AnandF11(c11[0],c11[1],jx,norm="none")
        self.F12 = AnandF12(c12[0],c12[1],jx,norm="none")
        self.F4F6 = AnandF4F6(c4[0],c4[1],c6[0],c6[1],jx,norm="none")
        self.F8F9 = AnandF8F9(c8[0],c8[1],c9[0],c9[1],jx,norm="none")
        self.F9F8 = AnandF9F8(c9[0],c9[1],c8[0],c8[1],jx,norm="none")
        self.F11F12 = AnandF11F12(c11[0],c11[1],c12[0],c12[1],jx,norm="none")
        self.F12F11 = AnandF12F11(c12[0],c12[1],c11[0],c11[1],jx,norm="none")
        
        self.FF = [self.F1,self.F4,self.F6,self.F7,self.F8,self.F9,self.F10,self.F11,self.F12,self.F4F6,self.F8F9,self.F9F8,self.F11F12,self.F12F11]

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        rate = 0
        for F in self.FF:
            rate+=F.dRdER(Target,ER,self.mX,self.sig,VelDist) # Sum the rate for each form factor
        return rate
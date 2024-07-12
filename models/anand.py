from dmmodel import DMModel
import numpy as np
from constants import *

"""
Class defintions for all the DM form factors, following the Anand formalism 
(i.e., unitless form factors as per arxiv 1308.6288)
These models are for use where you want to plot/constrain an experimental cross section
"""

class AnandF1(DMModel):
    def __init__(self, cp, cn):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self, Target, ER):
        """
        Form factor expression for this model
        """
        p_p = self.cp*self.cp*Target.FMpp(ER)
        p_n = self.cp*self.cn*Target.FMpn(ER)
        n_p = self.cn*self.cp*Target.FMnp(ER)
        n_n = self.cn*self.cn*Target.FMnn(ER)
        return p_p+p_n+n_p+n_n

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)


class AnandF3(DMModel):
    def __init__(self, cp, cn):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self, Target, ER,vm):
        """
        Form factor expression for this model
        """
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
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms)
            cross_sec = sig 
            dsigdER_g = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            dsigdER_h = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        

class AnandF4(DMModel):
    def __init__(self, cp, cn,jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self, Target, ER):
        """
        Form factor expression for this model
        """
        p_p = self.cp*self.cp*(Target.FS1pp(ER)+Target.FS2pp(ER))
        p_n = self.cp*self.cn*(Target.FS1pn(ER)+Target.FS2pn(ER))
        n_p = self.cn*self.cp*(Target.FS1np(ER)+Target.FS2np(ER))
        n_n = self.cn*self.cn*(Target.FS1nn(ER)+Target.FS2nn(ER))

        return Target.spin_dep(self.jx)*(p_p+p_n+n_p+n_n)/16
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target, ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        
class AnandF5(DMModel):
    def __init__(self, cp, cn, jx):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self, Target, ER, vm):
        """
        Form factor expression for this model
        """
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
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target, ER, vm/kms)
            cross_sec = sig 
            dsigdER_g = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            dsigdER_h = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        

class AnandF6(DMModel):
    def __init__(self, cp, cn, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self, Target, ER):
        """
        Form factor expression for this model
        """
        p_p = self.cp*self.cp*(Target.FS2pp(ER))
        p_n = self.cp*self.cn*(Target.FS2pn(ER))
        n_p = self.cn*self.cp*(Target.FS2np(ER))
        n_n = self.cn*self.cn*(Target.FS2nn(ER))

        return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,4.)*(p_p+p_n+n_p+n_n)/16

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class AnandF7(DMModel):
    def __init__(self, cp, cn):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self,Target,ER,vm):
        """
        O7,7 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O7 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        p_p = self.cp*self.cp*Target.FS1pp(ER)
        p_n = self.cp*self.cn*Target.FS1pn(ER)
        n_p = self.cn*self.cp*Target.FS1np(ER)
        n_n = self.cn*self.cn*Target.FS1nn(ER)

        h = (p_p+p_n+n_p+n_n)/8
        g = -vm*vm(p_p+p_n+n_p+n_n)/8

        return [g,h]

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms)
            cross_sec = sig 
            dsigdER_g = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            dsigdER_h = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        

class AnandF8(DMModel):
    def __init__(self, cp, cn, jx):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self,Target,ER,vm):
        """
        O7,7 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O7 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        g_p_p = self.cp*self.cp*(pow(Target.Q(ER)/mp,2)*Target.FDpp(ER)-vm*vm*Target.FMpp(ER))
        g_p_n = self.cp*self.cn*(pow(Target.Q(ER)/mp,2)*Target.FDpn(ER)-vm*vm*Target.FMpn(ER))
        g_n_p = self.cn*self.cp*(pow(Target.Q(ER)/mp,2)*Target.FDnp(ER)-vm*vm*Target.FMnp(ER))
        g_n_n = self.cn*self.cn*(pow(Target.Q(ER)/mp,2)*Target.FDnn(ER)-vm*vm*Target.FMnn(ER))
        g = 0.25*Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)
        
        h_p_p = self.cp*self.cp*Target.FMpp(ER)
        h_p_n = self.cp*self.cn*Target.FMpn(ER)
        h_n_p = self.cn*self.cp*Target.FMnp(ER)
        h_n_n = self.cn*self.cn*Target.FMnn(ER)
        h = 0.25*Target.spin_dep(self.jx)*(h_p_p+h_p_n+h_n_p+h_n_n)

        return [g,h]
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target, ER,vm/kms)
            cross_sec = sig 
            dsigdER_g = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            dsigdER_h = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        

class AnandF9(DMModel):
    def __init__(self, cp, cn, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self,Target,ER):
        """
        O7,7 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O7 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        p_p = self.cp*self.cp*Target.FS1pp(ER)
        p_n = self.cp*self.cn*Target.FS1pn(ER)
        n_p = self.cn*self.cp*Target.FS1np(ER)
        n_n = self.cn*self.cn*Target.FS1nn(ER)

        return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/16 

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class AnandF10(DMModel):
    def __init__(self, cp, cn):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self,Target,ER):
        """
        O7,7 operator, depends on DM spin
        cn and cp both couplings between relevant nucleon and O7 operator
        Annoyingly, this FF has both g and h depedence, so lets return a list
        Also depends on the min velocity of the DM
        """
        p_p = self.cp*self.cp*Target.FS2pp(ER)
        p_n = self.cp*self.cn*Target.FS2pn(ER)
        n_p = self.cn*self.cp*Target.FS2np(ER)
        n_n = self.cn*self.cn*Target.FS2nn(ER)

        return np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/4 

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class AnandF11(DMModel):
    def __init__(self, cp, cn,jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self, Target, ER):
        """
        Form factor expression for this model
        """
        p_p = self.cp*self.cp*Target.FMpp(ER)
        p_n = self.cp*self.cn*Target.FMpn(ER)
        n_p = self.cn*self.cp*Target.FMnp(ER)
        n_n = self.cn*self.cn*Target.FMnn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/4

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        
class AnandF12(DMModel):
    def __init__(self, cp, cn, jx):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

    def FF(self,Target,ER,vm): 
        """
        Defining O12,12 based on form factors
        cn and cp both couplings between relevant nucleon and O1 operator
        """
        
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

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms)
            cross_sec = sig 
            dsigdER_g = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            dsigdER_h = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        
###### Need to double check units and signs for these       
class AnandF1F3(DMModel):
    def __init__(self, c1p, c1n, c3p, c3n):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        c0 = np.sqrt(c1p**2 + c1n**2 + c3p**2 + c3n**2)
        self.c1p = c1p/c0
        self.c1n = c1n/c0
        self.c3p = c3p/c0
        self.c3n = c3n/c0

    def FF(self, Target, ER):
        """
        O1,3 operator, depends on DM spin
        """
        p_p = self.cp1*self.cp3*Target.FMPhi2pp(ER)
        p_n = self.cp1*self.cn3*Target.FMPhi2pn(ER)
        n_p = self.cn1*self.cp3*Target.FMPhi2np(ER)
        n_n = self.cn1*self.cn3*Target.FMPhi2nn(ER)
        return np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)   

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class AnandF4F5(DMModel):
    def __init__(self, c4p, c4n, c5p, c5n, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        c0 = np.sqrt(c4p**2 + c4n**2 + c5p**2 + c5n**2)
        self.c4p = c4p/c0
        self.c4n = c4n/c0
        self.c5p = c5p/c0
        self.c5n = c5n/c0
        self.jx = jx

    def FF(self, Target, ER):
        """
        O4,5 operator, depends on DM spin
        """
        p_p = self.cp4*self.cp5*Target.FS1Dpp(ER)
        p_n = self.cp4*self.cn5*Target.FS1Dpn(ER)
        n_p = self.cn4*self.cp5*Target.FS1Dnp(ER)
        n_n = self.cn4*self.cn5*Target.FS1Dnn(ER)
        return 0.25*Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class AnandF4F6(DMModel):
    def __init__(self, c4p, c4n, c6p, c6n, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        c0 = np.sqrt(c4p**2 + c4n**2 + c6p**2 + c6n**2)
        self.c4p = c4p/c0
        self.c4n = c4n/c0
        self.c6p = c6p/c0
        self.c6n = c6n/c0
        self.jx = jx

    def FF(self, Target, ER):
        """
        O4,5 operator, depends on DM spin
        """
        p_p = self.cp4*self.cp6*Target.FS2pp(ER)
        p_n = self.cp4*self.cn6*Target.FS2pn(ER)
        n_p = self.cn4*self.cp6*Target.FS2np(ER)
        n_n = self.cn4*self.cn6*Target.FS2nn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/mp,2.)*(p_p+p_n+n_p+n_n)/8

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            cross_sec = sig 
            dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
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
            FF = Target.F11(ER,self.cp,self.cn)
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
            FF = Target.F33(ER,self.cp,self.cn,vm/kms)
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
            FF = Target.F44(ER,self.cp,self.cn,self.jx)
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
            FF = Target.F55(ER,self.cp,self.cn,vm/kms,self.jx)
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
            FF = Target.F66(ER,self.cp,self.cn,self.jx)
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
            FF = Target.F77(ER,self.cp,self.cn,vm/kms)
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
            FF = Target.F88(ER,self.cp,self.cn,vm/kms,self.jx)
            cross_sec = sig 
            dsigdER_g = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            dsigdER_h = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        

class AnandF9(DMModel):
    def __init__(self, cp, cn,jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp = cp/np.sqrt(cp*cp+cn*cn)
        self.cn = cn/np.sqrt(cp*cp+cn*cn)
        self.jx = jx
        ## ultimately could try and use this to help with the mapping from EFT to exp
        ## eg, normalise them here to allow for a certain cross section 

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
            FF = Target.F99(ER,self.cp,self.cn,self.jx)
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
            FF = Target.F1010(ER,self.cp,self.cn)
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
            FF = Target.F1111(ER,self.cp,self.cn,self.jx)
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
            FF = Target.F1212(ER,self.cp,self.cn,vm/kms,self.jx)
            cross_sec = sig 
            dsigdER_g = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            dsigdER_h = cross_sec*FF[0]*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
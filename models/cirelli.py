from dmmodel import DMModel
import numpy as np
from constants import *

"""
Class defintions for all the DM form factors, following the Cirelli formalism 
(i.e., unitless form factors as per arxiv 1307.5955)
These models are for use where you want to plot/constrain a new physics scale
"""

######## NOTE TO SELF!
# Still need to compute a cross section properly for these
# Actually... the real high energy way to do it would be to pass the quark couplings.
# I.e., absorb the functions Raghda has worked on to compute the NR couplings given quark level...
# these functions are very much a work in progress.

class CirelliF11(DMModel):
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,lam,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [lam] = [eV] new physics scale

        Output units: cpd/kg/keV
        """
        vm = self.vmin(Target,mX,ER)
        cn_NR = 4*mX*mp*mp*0.433/pow(lam,3) # NR EFT coupling assuming EFT cn = 0.433/pow(lam*keV,3)
        cp_NR = 4*mX*mp*mp*0.3754/pow(lam,3) # NR EFT coupling assuming EFT cp = 0.3754/pow(lam*keV,3)
        FF = Target.F11(ER,cp_NR,cn_NR) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1. vmin needs to be unitless
        dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 ## units of 1/[eV]3, Cirelli cross section expression
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)

class CirelliF3(DMModel):
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
        

class CirelliF4(DMModel):
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
        
class CirelliF5(DMModel):
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
        

class CirelliF66(DMModel):
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,lam,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [lam] = [eV] new physics scale

        Output units: cpd/kg/keV
        """
        vm = self.vmin(Target,mX,ER)
        cp_NR = -4*1.91E-24
        cn_NR = 4*3.59E-25
        FF = pow(mp,4)*Target.F66(ER,cp_NR,cn_NR,0.5) ## need to add the mass term to get the unitful term
        
        dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 ## units of 1/[eV]3, Cirelli cross section expression
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class CirelliF7(DMModel):
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
        

class CirelliF8(DMModel):
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
        

class CirelliF9(DMModel):
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
        

class CirelliF10(DMModel):
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
        

class CirelliF11(DMModel):
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
        
class CirelliF12(DMModel):
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
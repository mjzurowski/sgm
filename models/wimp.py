from dmmodel import DMModel
import numpy as np
from constants import *

### Class definition for standard SI WIMP.
#### Note that this is the "blueprint" for DM models defined with NREFT. You should be able to just take this and switch out the FF terms

class SIWIMP(DMModel):
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
        FF = Target.F11(ER,1/np.sqrt(2),1/np.sqrt(2)) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = sig*2 ## as we normalise the nucleon vector [cn, cp] our built in cross section is sigma_N. Multiply by 2 to get this from sigma_p (which is generally assumed to be the input)
        dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
    
class WIMPO3(DMModel):
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
        vm = self.vmin(Target,mX,ER)
        FF = Target.F33(ER,1,1,vm/kms) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1. vmin needs to be unitless
        # maybe better to bake in the unit conversion on the form factor side?
        FF_g = FF[0] ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        FF_h = FF[1]
        cross_sec = sig ## as we normalise the nucleon vector [cn, cp] our built in cross section is sigma_N. Multiply by 2 to get this from sigma_p (which is generally assumed to be the input)
        dsigdER_g = cross_sec*FF_g*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        dsigdER_h = cross_sec*FF_h*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))


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
        FF = pow(mp,4)*Target.F66(ER,cp_NR,cn_NR,0.5) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1. vmin needs to be unitless
        
        dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 ## units of 1/[eV]3, Cirelli cross section expression
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
    
class AnandF11(DMModel):
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
        cp = 0.3754
        cn = 0.433
        FF = Target.F11(ER,cp/np.sqrt(cp*cp+cn*cn),cn/np.sqrt(cp*cp+cn*cn)) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = sig
        dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)


class AnandF66(DMModel):
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
        cp = -1.91E-24*mp/mX
        cn = 3.59E-25*mp/mX
        FF = Target.F66(ER,cp/np.sqrt(cp*cp+cn*cn),cn/np.sqrt(cp*cp+cn*cn),0.5) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = sig
        dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
######################################################
### Class definition for standard SI WIMP with Helm form factors
#### Nore that this is equivalent to NREFT with F11
    
class SIWIMP_Helm(DMModel):
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we take the Helm Form Factor rather than those defined from nrefts. Note that it should be ~ equal to the O1 SIWIMP
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: cpd/kg/keV
        """
        FF = Target.Helm(ER)**2 ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        dsigdER = (1/kg_to_eV)*sig*FF*Target.A()*Target.A()/(2*Target.N_T()*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
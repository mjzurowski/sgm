from dmmodel import DMModel
import numpy as np
from constants import *
import math
from Coefficient_Mapping import *

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
    
    def dRdER(self,Target,ER,mX,sig,VelDist,cp,cn):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        exposure = 612 #kg days
        #cpd_conversion_new = 2.597e-10
        cpd_conversion_new = 1014247.1936227841 # takes into account ev2Day# (1/cm2ev^3)
        FF = Target.F11(ER,cp,cn) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        cross_sec = sig*2 ## as we normalise the nucleon vector [cn, cp] our built in cross section is sigma_N. Multiply by 2 to get this from sigma_p (which is generally assumed to be the input)
        dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER)
        return exposure*cpd_conversion_new*Target.N_T()*(VelDist.rho*1e9/mX)*dsigdER*VelDist.gdist(vm)/((mp+mX)**2)


    def dRdER_New(self,Target,ER,mX,VelDist,cp,cn):
        exposure = 612 #kg days
        NA=6.022e26
        ev2Day=1.314e20
        kev2Day=1.314e23
        cm2eV=5.06e4
        ABar= 0.205*70 + 0.274*72 + 0.0776*73 + 0.365*74 + 0.0775*76
        FF = Target.F11(ER,cp,cn)
        dsigdER = FF*Target.mT() *NA/(32*math.pi*mX**3 *mp**2 *ABar)
        vm = self.vmin(Target,mX,ER)
        return exposure*kev2Day*(VelDist.rho*1e9/(cm2eV**3))*dsigdER*VelDist.gdist(vm)


class WIMPall(DMModel):
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)


    def dRdER_Old(self,Target,ER,mX,VelDist,Coeff,Lambda,jx):
        vm = self.vmin(Target,mX,ER)
        exposure = 612 #kg days
        NA=6.022e26
        ev2Day=1.314e20
        kev2Day=1.314e23
        cm2eV=5.06e4
        ABar= 0.205*70 + 0.274*72 + 0.0776*73 + 0.365*74 + 0.0775*76
        FF11 = Target.F11(ER,c1(p,Coeff,mX,Lambda,mp),c1(n,Coeff,mX,Lambda,mp))
        #FF33 = Target.F33(ER,c3(p,Coeff,mX,Lambda),c3(n,Coeff,mX,Lambda),vm/kms)
        FF44 = Target.F44(ER,c4(p,Coeff,mX,Lambda,mp),c4(n,Coeff,mX,Lambda,mp),jx)
        #FF55 =
        FF66 = Target.F66(ER,c6(p,Coeff,mX,Lambda,mp),c6(n,Coeff,mX,Lambda,mp),jx)
        FF77 = Target.F77(ER,c7(p,Coeff,mX,Lambda,mp),c7(n,Coeff,mX,Lambda,mp),vm/kms) #g,h
        FF88 = Target.F88(ER,c8(p,Coeff,mX,Lambda,mp),c8(n,Coeff,mX,Lambda,mp),jx,vm/kms) #g,h
        FF99 = Target.F99(ER,c9(p,Coeff,mX,Lambda,mp),c9(n,Coeff,mX,Lambda,mp),jx)
        FF1010 = Target.F1010(ER,c10(p,Coeff,mX,Lambda,mp),c10(n,Coeff,mX,Lambda,mp))
        FF1111 = Target.F1111(ER,c11(p,Coeff,mX,Lambda,mp),c11(n,Coeff,mX,Lambda,mp),jx)
        FF1212 = Target.F1212(ER,c12(p,Coeff,mX,Lambda,mp),c12(n,Coeff,mX,Lambda,mp),vm/kms,jx) #g,h
        #FF13
        #FF45
        FF46 = Target.F46(ER,c4(p,Coeff,mX,Lambda,mp),c4(n,Coeff,mX,Lambda,mp),c6(p,Coeff,mX,Lambda,mp),c6(n,Coeff,mX,Lambda,mp),jx)
        FF64 = Target.F64(ER,c6(p,Coeff,mX,Lambda,mp),c6(n,Coeff,mX,Lambda,mp),c4(p,Coeff,mX,Lambda,mp),c4(n,Coeff,mX,Lambda,mp),jx)
        FF98 = Target.F98(ER,c9(p,Coeff,mX,Lambda,mp),c9(n,Coeff,mX,Lambda,mp),c8(p,Coeff,mX,Lambda,mp),c8(n,Coeff,mX,Lambda,mp),jx)
        FF89 = Target.F89(ER,c8(p,Coeff,mX,Lambda,mp),c8(n,Coeff,mX,Lambda,mp),c9(p,Coeff,mX,Lambda,mp),c9(n,Coeff,mX,Lambda,mp),jx)
        FF1112 = Target.F1112(ER,c11(p,Coeff,mX,Lambda,mp),c11(n,Coeff,mX,Lambda,mp),c12(p,Coeff,mX,Lambda,mp),c12(n,Coeff,mX,Lambda,mp),jx)
        FF1211 = Target.F1211(ER,c12(p,Coeff,mX,Lambda,mp),c12(n,Coeff,mX,Lambda,mp),c11(p,Coeff,mX,Lambda,mp),c11(n,Coeff,mX,Lambda,mp),jx)
        dsigdER = Target.mT() *NA/(32*math.pi*mX**3 *mp**2 *ABar)
        return exposure*kev2Day*(VelDist.rho*1e9/(cm2eV**3))*dsigdER*((FF11+FF44+FF66+FF77[0]+FF88[0]+FF99+FF1010+FF1111+FF1212[0]+FF46+FF64+FF98+FF89+FF1112+FF1211)*VelDist.gdist(vm) + (FF77[1]+FF88[1]+FF1212[1])*VelDist.hdist(vm))

    
    def dRdER(self,Target,ER,mX,VelDist,Coeff,Lambda,jx):
        vm = self.vmin(Target,mX,ER)
        exposure = 612 #kg days
        NA=6.022e26
        ev2Day=1.314e20
        kev2Day=1.314e23
        cm2eV=5.06e4
        ABar= 0.205*70 + 0.274*72 + 0.0776*73 + 0.365*74 + 0.0775*76

        FF11=FF44=FF66=FF99=FF1010=FF1111=FF46=FF64=FF98=FF89=FF1112=FF1211=0
        FF77=FF88=FF1212=[0,0]
        
        if Coeff[0]==1:
            FF11 = Target.F11(ER,c1(p,Coeff,mX,Lambda,mp),c1(n,Coeff,mX,Lambda,mp))
        if Coeff[1]==1:
            FF1111 = Target.F1111(ER,c11(p,Coeff,mX,Lambda,mp),c11(n,Coeff,mX,Lambda,mp),jx)
        if Coeff[2]==1:
            FF1010 = Target.F1010(ER,c10(p,Coeff,mX,Lambda,mp),c10(n,Coeff,mX,Lambda,mp))
        if Coeff[3]==1:
            FF66 = Target.F66(ER,c6(p,Coeff,mX,Lambda,mp),c6(n,Coeff,mX,Lambda,mp),jx)
        if Coeff[4]==1:
            FF11 = Target.F11(ER,c1(p,Coeff,mX,Lambda,mp),c1(n,Coeff,mX,Lambda,mp))
        if Coeff[5]==1:
            FF88 = Target.F88(ER,c8(p,Coeff,mX,Lambda,mp),c8(n,Coeff,mX,Lambda,mp),jx,vm/kms) #g,h
            FF99 = Target.F99(ER,c9(p,Coeff,mX,Lambda,mp),c9(n,Coeff,mX,Lambda,mp),jx)
            FF98 = Target.F98(ER,c9(p,Coeff,mX,Lambda,mp),c9(n,Coeff,mX,Lambda,mp),c8(p,Coeff,mX,Lambda,mp),c8(n,Coeff,mX,Lambda,mp),jx)
            FF89 = Target.F89(ER,c8(p,Coeff,mX,Lambda,mp),c8(n,Coeff,mX,Lambda,mp),c9(p,Coeff,mX,Lambda,mp),c9(n,Coeff,mX,Lambda,mp),jx)
        if Coeff[6]==1:
            FF77 = Target.F77(ER,c7(p,Coeff,mX,Lambda,mp),c7(n,Coeff,mX,Lambda,mp),vm/kms) #g,h
            FF99 = Target.F99(ER,c9(p,Coeff,mX,Lambda,mp),c9(n,Coeff,mX,Lambda,mp),jx)
        if Coeff[7]==1:
            FF44 = Target.F44(ER,c4(p,Coeff,mX,Lambda,mp),c4(n,Coeff,mX,Lambda,mp),jx)
        if Coeff[8]==1:
            FF44 = Target.F44(ER,c4(p,Coeff,mX,Lambda,mp),c4(n,Coeff,mX,Lambda,mp),jx)
        if Coeff[9]==1:
            FF1010 = Target.F1010(ER,c10(p,Coeff,mX,Lambda,mp),c10(n,Coeff,mX,Lambda,mp))
            FF1111 = Target.F1111(ER,c11(p,Coeff,mX,Lambda,mp),c11(n,Coeff,mX,Lambda,mp),jx)
            FF1212 = Target.F1212(ER,c12(p,Coeff,mX,Lambda,mp),c12(n,Coeff,mX,Lambda,mp),vm/kms,jx) #g,h
            FF1112 = Target.F1112(ER,c11(p,Coeff,mX,Lambda,mp),c11(n,Coeff,mX,Lambda,mp),c12(p,Coeff,mX,Lambda,mp),c12(n,Coeff,mX,Lambda,mp),jx)
            FF1211 = Target.F1211(ER,c12(p,Coeff,mX,Lambda,mp),c12(n,Coeff,mX,Lambda,mp),c11(p,Coeff,mX,Lambda,mp),c11(n,Coeff,mX,Lambda,mp),jx)
        
        dsigdER = Target.mT() *NA/(32*math.pi*mX**3 *mp**2 *ABar)
        return exposure*kev2Day*(VelDist.rho*1e9/(cm2eV**3))*dsigdER*((FF11+FF44+FF66+FF77[0]+FF88[0]+FF99+FF1010+FF1111+FF1212[0]+FF46+FF64+FF98+FF89+FF1112+FF1211)*VelDist.gdist(vm) + (FF77[1]+FF88[1]+FF1212[1])*VelDist.hdist(vm))
      





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
       return (c*1.E-3)*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we take the Helm Form Factor rather than those defined from nrefts. Note that it should be ~ equal to the O1 SIWIMP
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        FF = Target.Helm(ER)**2 ## form factor with couplings. Note that proton and neutron couplings are normalised to 1
        dsigdER = (1/kg_to_eV)*sig*FF*Target.A()*Target.A()/(2*Target.N_T()*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
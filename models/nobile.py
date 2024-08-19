from dmmodel import DMModel
import numpy as np
from constants import *
from models.couplings import *

"""
Class defintions for all the DM form factors, following the Nobile formalism 
(i.e., unitful form factors as per arxiv 2104.12785)
These models are for use where you want to plot/constrain a new physics scale.
NB: not all of the NR operators correspond to a relativistic equivalent. We include them here for completeness, and allow the user to specify their own couplings.
"""

class NobileF1(DMModel):
    def __init__(self, cq):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c1 = c1_N(cq)
        self.c5 = c5_N(cq)

    def FF(self, Target, ER, mX, Lambda):
        """
        Form factor expression for O1
        Output units: unitless
        """
        cp = 4*mp*mX*(self.c1[0]/pow(Lambda,3)+self.c5[0]/pow(Lambda,2))
        cn = 4*mp*mX*(self.c1[1]/pow(Lambda,3)+self.c5[1]/pow(Lambda,2))
        p_p = cp*cp*Target.FMpp(ER)
        p_n = cp*cn*Target.FMpn(ER)
        n_p = cn*cp*Target.FMnp(ER)
        n_n = cn*cn*Target.FMnn(ER)
        return p_p+p_n+n_p+n_n

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c1==self.c5==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,mX,Lambda)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)

class NobileF3(DMModel):
    def __init__(self, cp, cn):
        """
        Initialise with a set of cp and cn values.
        The NR operator 3 has no high energy equivalent, so allow these to be defined manually
        """
        self.cp = cp
        self.cn = cn

    def FF(self, Target, ER,vm):
        """
        Form factor expression for O3
        """
        h_p_p = self.cp*self.cp*Target.FS1pp(ER)
        h_p_n = self.cp*self.cn*Target.FS1pn(ER)
        h_n_p = self.cn*self.cp*Target.FS1np(ER)
        h_n_n = self.cn*self.cn*Target.FS1nn(ER)
        h = np.power(Target.Q(ER),2)*(h_p_p+h_p_n+h_n_p+h_n_n)/8  #units = eV^2

        g_p_p = self.cp*self.cp*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2pp(ER)/4-np.power(vm*Target.Q(ER),2.)*Target.FS1pp(ER)/8)
        g_p_n = self.cp*self.cn*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2pn(ER)/4-np.power(vm*Target.Q(ER),2.)*Target.FS1pn(ER)/8)
        g_n_p = self.cn*self.cp*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2np(ER)/4-np.power(vm*Target.Q(ER),2.)*Target.FS1np(ER)/8)
        g_n_n = self.cn*self.cn*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2nn(ER)/4-np.power(vm*Target.Q(ER),2.)*Target.FS1nn(ER)/8)
        g = (g_p_p+g_p_n+g_n_p+g_n_n)  #units = eV^2

        return [g,h]
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        
class NobileF4(DMModel):
    def __init__(self, cq, jx):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c8 = c8_N(cq)
        self.c9 = c9_N(cq)
        self.jx = jx

    def FF(self, Target, ER, mX, Lambda):
        """
        Form factor expression for O4
        Output units: unitless
        """
        cp = 16*mp*mX*(2*self.c9[0]-self.c8[0])/pow(Lambda,2)
        cn = 16*mp*mX*(2*self.c9[1]-self.c8[1])/pow(Lambda,2)
        p_p = cp*cp*(Target.FS1pp(ER)+Target.FS2pp(ER))
        p_n = cp*cn*(Target.FS1pn(ER)+Target.FS2pn(ER))
        n_p = cn*cp*(Target.FS1np(ER)+Target.FS2np(ER))
        n_n = cn*cn*(Target.FS1nn(ER)+Target.FS2nn(ER))

        return Target.spin_dep(self.jx)*(p_p+p_n+n_p+n_n)/16
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c8==self.c9==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target, ER, mX, Lambda)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        
class NobileF5(DMModel):
    def __init__(self, cp, cn, jx):
        """
        Initialise with a set of cp and cn values
        The NR operator 3 has no high energy equivalent, so allow these to be defined manually
        """
        self.cp = cp
        self.cn = cn
        self.jx = jx

    def FF(self, Target, ER, vm):
        """
        Form factor expression for O5
        """
        h_p_p = self.cp*self.cp*Target.FMpp(ER)
        h_p_n = self.cp*self.cn*Target.FMpn(ER)
        h_n_p = self.cn*self.cp*Target.FMnp(ER)
        h_n_n = self.cn*self.cn*Target.FMnn(ER)
        h = Target.spin_dep(self.jx)*np.power(Target.Q(ER),2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/4 #units = eV^2

        g_p_p = self.cp*self.cp*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FDpp(ER)-np.power(vm*Target.Q(ER),2.)*Target.FMpp(ER))
        g_p_n = self.cp*self.cn*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FDpn(ER)-np.power(vm*Target.Q(ER),2.)*Target.FMpn(ER))
        g_n_p = self.cn*self.cp*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FDnp(ER)-np.power(vm*Target.Q(ER),2.)*Target.FMnp(ER))
        g_n_n = self.cn*self.cn*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FDnn(ER)-np.power(vm*Target.Q(ER),2.)*Target.FMnn(ER))
        g = Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/4 #units = eV^2

        return [g,h]
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target, ER, vm/kms)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        
class NobileF6(DMModel):
    def __init__(self, cq, jx):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c4 = c4_N(cq)
        self.jx = jx

    def FF(self, Target, ER, Lambda):
        """
        Form factor expression for O6

        Output units: unitless
        """
        cp = 4*self.c4[0]/pow(Lambda,3)
        cn = 4*self.c4[1]/pow(Lambda,3)
        p_p = cp*cp*(Target.FS2pp(ER))
        p_n = cp*cn*(Target.FS2pn(ER))
        n_p = cn*cp*(Target.FS2np(ER))
        n_n = cn*cn*(Target.FS2nn(ER))

        return Target.spin_dep(self.jx)*np.power(Target.Q(ER),4.)*(p_p+p_n+n_p+n_n)/16

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: [cm^2]/[eV] 
        """
        if(self.c4==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,Lambda)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        
class NobileF7(DMModel):
    def __init__(self, cq):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c7 = c7_N(cq)

    def FF(self,Target,ER,vm,mX,Lambda):
        """
        Form factor expression for O7
        """
        cp = -8*mp*mX*self.c7[0]/pow(Lambda,2)
        cn = -8*mp*mX*self.c7[1]/pow(Lambda,2)
        p_p = cp*cp*Target.FS1pp(ER)
        p_n = cp*cn*Target.FS1pn(ER)
        n_p = cn*cp*Target.FS1np(ER)
        n_n = cn*cn*Target.FS1nn(ER)

        h = (p_p+p_n+n_p+n_n)/8
        g = -vm*vm*(p_p+p_n+n_p+n_n)/8

        return [g,h]

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c7==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms,mX,Lambda)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        
class NobileF8(DMModel):
    def __init__(self, cq, jx):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c6 = c6_N(cq)
        self.jx = jx

    def FF(self,Target,ER,vm,mX,Lambda):
        """
        Form factor expression for O8

        Output units: unitless
        """
        cp = 8*mp*mX*self.c6[0]/pow(Lambda,2)
        cn = 8*mp*mX*self.c6[1]/pow(Lambda,2)

        g_p_p = cp*cp*(pow(Target.Q(ER)/mp,2)*Target.FDpp(ER)-vm*vm*Target.FMpp(ER))
        g_p_n = cp*cn*(pow(Target.Q(ER)/mp,2)*Target.FDpn(ER)-vm*vm*Target.FMpn(ER))
        g_n_p = cn*cp*(pow(Target.Q(ER)/mp,2)*Target.FDnp(ER)-vm*vm*Target.FMnp(ER))
        g_n_n = cn*cn*(pow(Target.Q(ER)/mp,2)*Target.FDnn(ER)-vm*vm*Target.FMnn(ER))
        g = 0.25*Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)
        
        h_p_p = cp*cp*Target.FMpp(ER)
        h_p_n = cp*cn*Target.FMpn(ER)
        h_n_p = cn*cp*Target.FMnp(ER)
        h_n_n = cn*cn*Target.FMnn(ER)
        h = 0.25*Target.spin_dep(self.jx)*(h_p_p+h_p_n+h_n_p+h_n_n)

        return [g,h]
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c6==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target, ER,vm/kms,mX,Lambda)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
        
class NobileF9(DMModel):
    def __init__(self, cq, jx):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c6 = c6_N(cq)
        self.c7 = c7_N(cq)
        self.jx = jx

    def FF(self,Target,ER,mX,Lambda):
        """
        Form factor expression for O9

        Output units: unitless
        """
        cp = 8*mX*self.c6[0]/pow(Lambda,2)+8*mp*self.c7[0]/pow(Lambda,2)
        cn = 8*mX*self.c6[1]/pow(Lambda,2)+8*mp*self.c7[1]/pow(Lambda,2)
        p_p = cp*cp*Target.FS1pp(ER)
        p_n = cp*cn*Target.FS1pn(ER)
        n_p = cn*cp*Target.FS1np(ER)
        n_n = cn*cn*Target.FS1nn(ER)

        return Target.spin_dep(self.jx)*np.power(Target.Q(ER),2.)*(p_p+p_n+n_p+n_n)/16 

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c6==self.c7==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,mX,Lambda)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        
class NobileF10(DMModel):
    def __init__(self, cq):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c3 = c3_N(cq)
        self.c10 = c10_N(cq)

    def FF(self,Target,ER,mX,Lambda):
        """
        Form factor expression for O10

        Output units: unitless
        """
        cp = 4*mX*self.c3[0]/pow(Lambda,3)-8*mp*self.c10[0]/pow(Lambda,2)
        cn = 4*mX*self.c3[1]/pow(Lambda,3)-8*mp*self.c10[1]/pow(Lambda,2)

        p_p = cp*cp*Target.FS2pp(ER)
        p_n = cp*cn*Target.FS2pn(ER)
        n_p = cn*cp*Target.FS2np(ER)
        n_n = cn*cn*Target.FS2nn(ER)

        return np.power(Target.Q(ER),2.)*(p_p+p_n+n_p+n_n)/4 

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c3==self.c10==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,mX,Lambda)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        
class NobileF11(DMModel):
    def __init__(self, cq, jx):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c2 = c2_N(cq)
        self.c10 = c10_N(cq)
        self.jx = jx

    def FF(self, Target, ER, mX, Lambda):
        """
        Form factor expression for O11

        Output units: unitless
        """
        cp = -4*mp*self.c2[0]/pow(Lambda,3)+8*mX*self.c10[0]/pow(Lambda,2)
        cn = -4*mp*self.c2[1]/pow(Lambda,3)+8*mX*self.c10[1]/pow(Lambda,2)

        p_p = cp*cp*Target.FMpp(ER)
        p_n = cp*cn*Target.FMpn(ER)
        n_p = cn*cp*Target.FMnp(ER)
        n_n = cn*cn*Target.FMnn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER),2.)*(p_p+p_n+n_p+n_n)/4 

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c2==self.c10==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,mX,Lambda)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        
class NobileF12(DMModel):
    def __init__(self, cq, jx):
        """
        Initialise with a set of quark couplings and define the relativistic ones needed for this model
        """
        self.c10 = c10_N(cq)
        self.jx = jx

    def FF(self,Target,ER,vm,mX,Lambda): 
        """
        Form factor expression for O12

        Output units: unitless
        """
        cp = -32*mp*mX*self.c10[0]/pow(Lambda,2)
        cn = -32*mp*mX*self.c10[1]/pow(Lambda,2)

        h_p_p = cp*cp*(Target.FS1pp(ER)/2 + Target.FS2pp(ER))
        h_p_n = cp*cn*(Target.FS1pn(ER)/2 + Target.FS2pn(ER))
        h_n_p = cn*cp*(Target.FS1np(ER)/2 + Target.FS2np(ER))
        h_n_n = cn*cn*(Target.FS1nn(ER)/2 + Target.FS2nn(ER))
        h = Target.spin_dep(self.jx)*(h_p_p+h_p_n+h_n_p+h_n_n)/16 

        g_p_p = cp*cp*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2pp(ER) + Target.FPhipp(ER))-np.power(vm,2.)*(Target.FS1pp(ER)/2 + Target.FS2pp(ER)))
        g_p_n = cp*cn*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2pn(ER) + Target.FPhipn(ER))-np.power(vm,2.)*(Target.FS1pn(ER)/2 + Target.FS2pn(ER)))
        g_n_p = cn*cp*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2np(ER) + Target.FPhinp(ER))-np.power(vm,2.)*(Target.FS1np(ER)/2 + Target.FS2np(ER)))
        g_n_n = cn*cn*(np.power(Target.Q(ER)/mp,2.)*(Target.FPhi2nn(ER) + Target.FPhinn(ER))-np.power(vm,2.)*(Target.FS1nn(ER)/2 + Target.FS2nn(ER)))
        g = Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  

        return [g,h]

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,Lambda,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.c10==[0,0]):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms,mX,Lambda)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2 # units of [cm]^2/[eV]
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))

class NobileF13(DMModel):
    def __init__(self, cp, cn, jx):
        """
        Initialise with a set of cp and cn values
        """
        self.cp = cp
        self.cn = cn
        self.jx = jx

    def FF(self,Target,ER,vm): 
        """
        Form factor expression for O13
        """
        h_p_p = self.cp*self.cp*(Target.FS2pp(ER))
        h_p_n = self.cp*self.cn*(Target.FS2pn(ER))
        h_n_p = self.cn*self.cp*(Target.FS2np(ER))
        h_n_n = self.cn*self.cn*(Target.FS2nn(ER))
        h = Target.spin_dep(self.jx)*np.power(Target.Q(ER),2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/16 # units = eV^2

        g_p_p = self.cp*self.cp*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhipp(ER)-np.power(vm*Target.Q(ER),2.)*self.FS2pp(ER))
        g_p_n = self.cp*self.cn*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhipn(ER)-np.power(vm*Target.Q(ER),2.)*self.FS2pn(ER))
        g_n_p = self.cn*self.cp*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhinp(ER)-np.power(vm*Target.Q(ER),2.)*self.FS2np(ER))
        g_n_n = self.cn*self.cn*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhinn(ER)-np.power(vm*Target.Q(ER),2.)*self.FS2nn(ER))
        g = Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  #units = eV^2

        return [g,h]

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))

class NobileF14(DMModel):
    def __init__(self, cp, cn):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp = cp
        self.cn = cn

    def FF(self, Target, ER,vm):
        """
        Form factor expression for O14
        """
        p_p = self.cp*self.cp*(Target.FS1pp(ER))
        p_n = self.cp*self.cn*(Target.FS1pn(ER))
        n_p = self.cn*self.cp*(Target.FS1np(ER))
        n_n = self.cn*self.cn*(Target.FS1nn(ER))

        h = np.power(Target.Q(ER),2.)*(p_p+p_n+n_p+n_n)/32    # units = eV^2
        g = -vm*vm*np.power(Target.Q(ER),2.)*(p_p+p_n+n_p+n_n)/32   # units = eV^2

        return [g,h]
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
 
class NobileF15(DMModel):
    def __init__(self, cp, cn):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp = cp
        self.cn = cn

    def FF(self, Target, ER,vm):
        """
        Form factor expression for O15
        """
        h_p_p = self.cp*self.cp*(Target.FS1pp(ER))
        h_p_n = self.cp*self.cn*(Target.FS1pn(ER))
        h_n_p = self.cn*self.cp*(Target.FS1np(ER))
        h_n_n = self.cn*self.cn*(Target.FS1nn(ER))

        h = Target.spin_dep(self.jx)*np.power(Target.Q(ER),4.)*(h_p_p+h_p_n+h_n_p+h_n_n)/32    # units = eV^4

        g_p_p = self.cp*self.cp*(np.power(Target.Q(ER)/(mp**(1/3)),6.)*Target.FPhi2pp(ER)-np.power(np.sqrt(vm)*Target.Q(ER),4.)*self.FS1pp(ER)/2)
        g_p_n = self.cp*self.cn*(np.power(Target.Q(ER)/(mp**(1/3)),6.)*Target.FPhi2pn(ER)-np.power(np.sqrt(vm)*Target.Q(ER),4.)*self.FS1pn(ER)/2)
        g_n_p = self.cn*self.cp*(np.power(Target.Q(ER)/(mp**(1/3)),6.)*Target.FPhi2np(ER)-np.power(np.sqrt(vm)*Target.Q(ER),4.)*self.FS1np(ER)/2)
        g_n_n = self.cn*self.cn*(np.power(Target.Q(ER)/(mp**(1/3)),6.)*Target.FPhi2nn(ER)-np.power(np.sqrt(vm)*Target.Q(ER),4.)*self.FS1nn(ER)/2)
        g = Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  #units=ev^4

        return [g,h]
    
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER,vm/kms)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
 

        
###### Need to double check units and signs for these       
class NobileF1F3(DMModel):
    def __init__(self, cp1, cn1, cp3, cn3):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp1 = cp1
        self.cn1 = cn1
        self.cp3 = cp3
        self.cn3 = cn3

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O1 and O3
        """
        p_p = self.cp1*self.cp3*Target.FMPhi2pp(ER)
        p_n = self.cp1*self.cn3*Target.FMPhi2pn(ER)
        n_p = self.cn1*self.cp3*Target.FMPhi2np(ER)
        n_n = self.cn1*self.cn3*Target.FMPhi2nn(ER)
        return -0.5*np.power(Target.Q(ER)/np.sqrt(mp),2.)*(p_p+p_n+n_p+n_n) # units = eV

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF4F5(DMModel):
    def __init__(self, cp4, cn4, cp5, cn5, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp4 = cp4
        self.cn4 = cn4
        self.cp5 = cp5
        self.cn5 = cn5
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O4 and O5
        """
        p_p = self.cp4*self.cp5*Target.FS1Dpp(ER)
        p_n = self.cp4*self.cn5*Target.FS1Dpn(ER)
        n_p = self.cn4*self.cp5*Target.FS1Dnp(ER)
        n_n = self.cn4*self.cn5*Target.FS1Dnn(ER)
        return -Target.spin_dep(self.jx)*np.power(Target.Q(ER)/np.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF4F6(DMModel):
    def __init__(self, cp4, cn4, cp6, cn6, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp4 = cp4
        self.cn4 = cn4
        self.cp6 = cp6
        self.cn6 = cn6
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O4 and O5
        """
        p_p = self.cp4*self.cp6*Target.FS2pp(ER)
        p_n = self.cp4*self.cn6*Target.FS2pn(ER)
        n_p = self.cn4*self.cp6*Target.FS2np(ER)
        n_n = self.cn4*self.cn6*Target.FS2nn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER),2.)*(p_p+p_n+n_p+n_n)/16 # units = eV^2

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF6F4(DMModel):
    def __init__(self, cp4, cn4, cp6, cn6, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp4 = cp4
        self.cn4 = cn4
        self.cp6 = cp6
        self.cn6 = cn6
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O4 and O6
        """
        p_p = self.cp6*self.cp4*Target.FS2pp(ER)
        p_n = self.cp6*self.cn4*Target.FS2pn(ER)
        n_p = self.cn6*self.cp4*Target.FS2np(ER)
        n_n = self.cn6*self.cn4*Target.FS2nn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER),2.)*(p_p+p_n+n_p+n_n)/16 # units = eV^2

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF8F9(DMModel):
    def __init__(self, cp8, cn8, cp9, cn9, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp8 = cp8
        self.cn8 = cn8
        self.cp9 = cp9
        self.cn9 = cn9
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O8 and O9
        """
        p_p = self.cp8*self.cp9*Target.FS1Dpp(ER)
        p_n = self.cp8*self.cn9*Target.FS1Dpn(ER)
        n_p = self.cn8*self.cp9*Target.FS1Dnp(ER)
        n_n = self.cn8*self.cn9*Target.FS1Dnn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/np.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF9F8(DMModel):
    def __init__(self, cp8, cn8, cp9, cn9, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp8 = cp8
        self.cn8 = cn8
        self.cp9 = cp9
        self.cn9 = cn9
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O8 and O9
        """
        p_p = self.cp9*self.cp8*Target.FS1Dpp(ER)
        p_n = self.cp9*self.cn8*Target.FS1Dpn(ER)
        n_p = self.cn9*self.cp8*Target.FS1Dnp(ER)
        n_n = self.cn9*self.cn8*Target.FS1Dnn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER)/np.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF11F12(DMModel):
    def __init__(self, cp11, cn11, cp12, cn12, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp11 = cp11
        self.cn11 = cn11
        self.cp12 = cp12
        self.cn12 = cn12
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O11 and O12
        """
        p_p = self.cp11*self.cp12*Target.FMPhi2pp(ER)
        p_n = self.cp11*self.cn12*Target.FMPhi2pn(ER)
        n_p = self.cn11*self.cp12*Target.FMPhi2np(ER)
        n_n = self.cn11*self.cn12*Target.FMPhi2nn(ER)
        return -Target.spin_dep(self.jx)*np.power(Target.Q(ER)/np.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF12F11(DMModel):
    def __init__(self, cp11, cn11, cp12, cn12, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp11 = cp11
        self.cn11 = cn11
        self.cp12 = cp12
        self.cn12 = cn12
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O12 and O11
        """
        p_p = self.cp12*self.cp11*Target.FMPhi2pp(ER)
        p_n = self.cp12*self.cn11*Target.FMPhi2pn(ER)
        n_p = self.cn12*self.cp11*Target.FMPhi2np(ER)
        n_n = self.cn12*self.cn11*Target.FMPhi2nn(ER)
        return -Target.spin_dep(self.jx)*np.power(Target.Q(ER)/np.sqrt(mp),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF11F15(DMModel):
    def __init__(self, cp11, cn11, cp15, cn15, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp11 = cp11
        self.cn11 = cn11
        self.cp15 = cp15
        self.cn15 = cn15
        self.jx = jx

    def FF(self, Target, ER):
        """
        Form factor expression for interference of O11 and O15
        """
        p_p = self.cp11*self.cp15*Target.FMPhi2pp(ER)
        p_n = self.cp11*self.cn15*Target.FMPhi2pn(ER)
        n_p = self.cn11*self.cp15*Target.FMPhi2np(ER)
        n_n = self.cn11*self.cn15*Target.FMPhi2nn(ER)
        return Target.spin_dep(self.jx)*np.power(Target.Q(ER)**2/(mp**(1/2)),2.)*(p_p+p_n+n_p+n_n)/8 # units = eV^3

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER = FF*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)
        

class NobileF12F15(DMModel):
    def __init__(self, cp12, cn12, cp15, cn15, jx):
        """
        Initialise with a set of cp and cn values (coupling to n and p) and DM spin
        """
        self.cp12 = cp12
        self.cn12 = cn12
        self.cp15 = cp15
        self.cn15 = cn15
        self.jx = jx

    def FF(self, Target, ER, vm):
        """
        Form factor expression for interference of O12 and O15
        """
        h_p_p = self.cp12*self.cp15*Target.FS1pp(ER)
        h_p_n = self.cp12*self.cn15*Target.FS1pn(ER)
        h_n_p = self.cn12*self.cp15*Target.FS1np(ER)
        h_n_n = self.cn12*self.cn15*Target.FS1nn(ER)

        h = -Target.spin_dep(self.jx)*np.power(Target.Q(ER),2.)*(h_p_p+h_p_n+h_n_p+h_n_n)/32  #units = eV^2

        g_p_p = self.cp12*self.cp15*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2pp(ER)-np.power(vm*Target.Q(ER),2.)*self.FS1pp(ER)/2)
        g_p_n = self.cp12*self.cp15*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2pn(ER)-np.power(vm*Target.Q(ER),2.)*self.FS1pn(ER)/2)
        g_n_p = self.cn12*self.cp15*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2np(ER)-np.power(vm*Target.Q(ER),2.)*self.FS1np(ER)/2)
        g_n_n = self.cn12*self.cn15*(np.power(Target.Q(ER)/np.sqrt(mp),4.)*Target.FPhi2nn(ER)-np.power(vm*Target.Q(ER),2.)*self.FS1nn(ER)/2)
        g = -Target.spin_dep(self.jx)*(g_p_p+g_p_n+g_n_p+g_n_n)/16  #units = eV^2

        return [g,h]    

    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return kms*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy

        Output units: cpd/kg/keV
        """
        if(self.cn==self.cp==0):
            # both coupling constants are zero, so the rate will be too
            return 0
        else:
            vm = self.vmin(Target,mX,ER)
            FF = self.FF(Target,ER)
            dsigdER_g = FF[0]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            dsigdER_h = FF[1]*Target.mT()/(32*np.pi*mX*mX*mp*mp)*eV2_to_cm2
            return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(dsigdER_g*VelDist.gdist(vm) + dsigdER_h*VelDist.hdist(vm))
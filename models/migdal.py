from dmmodel import DMModel
import numpy as np
import scipy.integrate as integrate
from constants import *

### Class definition for MIGDAL Effect using Standard WIMP Assumptions.

class MIGDAL(DMModel):

    def vmin(self, Target, mX, ER, EM):
        
        """
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [EE] = [eV] Ionization energy (Electron + Binding)

        Output units: [km/s]

        This should return an array of vmin values accounting for each of the transition elements.x
        """
        
        vmins = [] # VMIN is returned as a list based on the transition energies in the target.
        for Eb in Target.eNL():
            vmins.append(c*1e-3*(Target.mT()*ER+Target.mu_T(mX)*(EM-Eb))/(Target.mu_T(mX)*np.power(2*Target.mT()*ER,0.5)))
        
        return vmins
        # return c*1e-3*(Target.mT()*ER+Target.mu_T(mX)*(EM))/(Target.mu_T(mX)*np.power(2*Target.mT()*ER,0.5))
 
    def dR2dEREE(self, Target, mX, ER, sig, VelDist, EM, SI = True, cp=1/np.sqrt(2), cn = 1/np.sqrt(2), jx = 1/2, Helm = True):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX]  = [eV] DM mass
        [ER]  = [eV] DM recoil energy
        [EM]  = [eV] Ionization energy (Electron + Binding)
        [sig] = [cm]^2 cross section
        [EE]  = [eV] Ionization energy (Electron + Binding)

        Output units: [cm^2]/[eV] 
        """

        if SI:
            if Helm:
                FF = Target.Helm(ER)**2
            else:
                FF = Target.F11(ER, cn=1/np.sqrt(2), cp=1/np.sqrt(2))
        else:
            FF = Target.F44(ER, cn,cp, jx)

        dsigdER = (1/kg_to_eV)*sig*FF*Target.A()**2 / (2*Target.N_T()*Target.mu_T(mX)**2) # units of [cm^2]/[eV]
        prefac = cpd_conversion*Target.N_T()*dsigdER*VelDist.rho/mX                       # 

        vmins = self.vmin(Target, mX, ER, EM)   # [km/s]       
        gdists = VelDist.gdist(vmins)           # Similarly taking and evaluating for the relevant vmin. [unitless]
        probs  = np.array(Target.eTrans((EM)))  # [unitless]

        # return prefac*np.array(gdists*probs)
        return np.sum(prefac*np.array(gdists*probs))
        # return (prefac*np.array(gdists*probs), vmins, gdists, probs)


    def dRdER(self, Target, E, mX, sig, VelDist, ERlim = 2*keV, **kwargs):
        return integrate.quad(lambda ER: self.dR2dEREE(Target, 
                                                        mX, 
                                                        ER, 
                                                        sig, 
                                                        VelDist,
                                                        EM=E, 
                                                        **kwargs), 0, ERlim)[0]

class SD_MIGDAL(DMModel): 
    # This is a test branch/sanity check with this implementation we dont account for/need to include the Eem energy with that of the vmin. 
    def vmin(self,Target,mX,ER):
       """
       [mX] = [eV] DM mass
       [ER] = [eV] DM recoil energy

       Output units: [km/s]
       """
       return (c*1.E-3)*np.abs((Target.mT()*ER/Target.mu_T(mX)))/np.power(2.*Target.mT()*ER,0.5)
    
    def dRdER(self,Target,ER,mX,sig,VelDist):
        """
        For this model, we just take coupling of n and p to be equal, and the only operator we care about is O1
        [mX] = [eV] DM mass
        [ER] = [eV] DM recoil energy
        [sig] = [cm]^2 cross section

        Output units: [cm^2]/[eV] 
        """
        # FF = Target.F11(ER,1/np.sqrt(2),1/np.sqrt(2)) ## form factor with couplings. Note that proton and neutron couplings are normalised to 1\
        FF = Target.F44(ER, 1, 0, 0.5)
        cross_sec = sig*2 ## as we normalise the nucleon vector [cn, cp] our built in cross section is sigma_N. Multiply by 2 to get this from sigma_p (which is generally assumed to be the input)
        dsigdER = cross_sec*FF*Target.mT()/(2*Target.mu_N(mX)*Target.mu_N(mX)) ## units of [cm^2]/[eV]
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*dsigdER*VelDist.gdist(vm)

    def dRdE(self, Target, mX, sig, VelDist, EE):
        return integrate.quad(lambda ER: self.dRdER(Target, ER,mX, sig, VelDist), a=0, b=np.inf)[0]*(1/np.pi)*np.array(Target.eTrans(EE))
    
from dmmodel import DMModel
import numpy as np
import scipy.integrate as integrate
from constants import *


class Migdal(DMModel):
    def vmin(self, Target, mX, ENR, EEM):
        vmin = c * 1e-3 * (np.sqrt(Target.mT() * ENR / (2. * Target.mu_T(mX)**2)) + EEM / np.sqrt(2 * Target.mT() * ENR))
        
        return vmin
 
    def dR2dENRdER(self, ENR, ER, Target, mX, sig, VelDist, SI=True, cp=1/np.sqrt(2), cn=1/np.sqrt(2), jx=1/2):
        if SI:
            FF = Target.F11(ENR, cn=1/np.sqrt(2), cp=1/np.sqrt(2))
        else:
            FF = Target.F44(ENR, cn,cp, jx)

        dsigmadENR = 2. * sig * FF * Target.mT() / (2 * Target.mu_N(mX)**2)
        prefactor = cpd_conversion * Target.N_T() * VelDist.rho / mX

        if Target.A() == 23:
            EEM = ER - (ENR * Target.quenching_factor())
        else:
            EEM = ER - (ENR * Target.quenching_factor())

        vmin = self.vmin(Target, mX, ENR, EEM) 
        gdist = VelDist.gdist(vmin)

        scaling =  2 * me**2 * ENR / (Target.mT())
        probs  = np.array(Target.eTrans((EEM))) / (2 * np.pi) * scaling

        return np.sum(prefactor * dsigmadENR * np.array(gdist * probs))

    def dRdER(self, Target, ER, mX, sig, VelDist, **kwargs):
        vmax_km_s = VelDist.vesc_km_s() + earth_v_km_s
        ER_max_eV = 2. * (vmax_km_s / (c * 1e-3) * Target.mu_T(mX))**2 / Target.mT()
        
        energy_arr = np.linspace(0.01, ER_max_eV, 100)
        rate_arr = [self.dR2dENRdER(ENR, ER, Target, mX, sig, VelDist, **kwargs) for ENR in energy_arr]

        return integrate.trapz(rate_arr, energy_arr)
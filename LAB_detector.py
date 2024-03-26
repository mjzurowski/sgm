from abc import ABC, abstractmethod
import scipy as sp
import numpy as np
import scipy.integrate as integrate
from constants import *

def quenching_proton(T):
        weight_H = 1 * (2 * 10+6) / (12 * (6 + 10) + 1 * (2 * 10 + 6))
        weight_C12 = 12 * (6 + 10) / (12 * (6 + 10) + 1 * (2 * 10 + 6))
        density_LAB = 0.859 # [g/cm3]

        eloss_file_H=np.loadtxt('models/neutrinos/p_in_H.dat', skiprows=8)
        eloss_file_C12=np.loadtxt('models/neutrinos/p_in_C12.dat', skiprows=8)
        eloss_file_energy = eloss_file_H[:, 0]
        eloss_file_LAB = (weight_H * eloss_file_H[:,1] + weight_C12 * eloss_file_C12[:,1]) * density_LAB

        integration_range_edges = np.geomspace(1e-4, T, 1000)
        integration_range_centers = 0.5 * (integration_range_edges[1:] + integration_range_edges[:-1])
        integration_range_diffs = np.diff(integration_range_edges)
        dEdx_interp = np.interp(integration_range_centers, eloss_file_energy, eloss_file_LAB, right=0)
        
        integral = np.sum([1 / (1 + 0.01 * dEdx) * dE for dEdx, dE in zip(dEdx_interp, integration_range_diffs)])
        
        return integral

class LABDetector(ABC):
    recoil_Es = np.geomspace(0.01, 100, 100)
    observed_Es = [quenching_proton(E) for E in recoil_Es]

    def ER_E(self, E):
        """
        Recoil energy as a function of observed energy
        (e.g., quenching factor, ionisation yield)
        """
        ER_E_interp = np.interp(E, self.observed_Es, self.recoil_Es)

        return ER_E_interp

    def dERdE(self, E):
        """
        Derivative of recoil energy wrt observed energy
        """
        dERdE = np.gradient(self.recoil_Es) / np.gradient(self.observed_Es)
        dERdE_interp = np.interp(E, self.observed_Es, dERdE)

        return dERdE_interp

    def dRdE_True(self, E, Func, NR=False,
                  per_proton=False, per_electron=False, per_carbon=False,
                  **kwargs):
        """
        Rate as a function of observed energy
        Assume "Func" object is something like a background or signal model that is a function of E and target. (NB background will probably be independent of target so it will be a useless variable there)
        """
        ratio_HC=1.639 # https://arxiv.org/pdf/1507.05613.pdf page163
        Np = 1e6 * ratio_HC * 1.007 / (ratio_HC * 1.007 + 12.011) / 1.007 * 6.02e23
        NC = 1e6 * 12.011 / (ratio_HC * 1.007 + 12.011) / 12.011 * 6.02e23
        Ne = Np + 6 * NC

        number_scaling = None
        if per_proton:
            number_scaling = Np
        elif per_electron:
            number_scaling = Ne
        elif per_carbon:
            number_scaling = NC

        if NR:
            ER = self.ER_E(E / MeV) * MeV
            dERdE = self.dERdE(E / MeV)
            Rate = dERdE * float(Func(ER, **kwargs)) * number_scaling
        else:
            Rate = float(Func(E, **kwargs)) * number_scaling

        return Rate

### General detector terms
    @abstractmethod
    def Emax(self):
        """
        Max energy of operation (used to define range of integrals)
        """
        pass

    @abstractmethod
    def Eff(self,E):
        """
        Efficiency as a function of observed energy
        """
        pass

    @abstractmethod
    def DeltaE(self,E):
        """
        Resolution as a function of observed energy
        """
        pass

    @abstractmethod
    def Res(self,E1,E2):
        """
        Resolution smearing function (probably a gaussian but lets leave arbitrary for now...)
        We assume E1 is the observed energy (E' in accompanying documentation) and E2 is the energy that will be integrated over (E_ee in accompanying documentation)
        """
        pass

    def dRdE(self, E, Model, NR=True, DE=0.01, Emax=None, **kwargs):
        """
        Observed rate for Model object smeared with resolution
        NR is a flag (default set true) that allows you to turn on and off nuclear recoil vs electron recoil (which will have different conversion factors to observed energy)
        DE is the step size for the energy array used for integration
        """
        if Emax is None:
            Emax = 2. * self.Emax()
        energy_arr = np.arange(0.01, Emax, DE)
        rate_arr = [self.dRdE_True(E2,Model,NR,**kwargs)*self.Res(E,E2) for E2 in energy_arr]
        return integrate.trapz(rate_arr,energy_arr)

    
import numpy as np
from constants import *

#### Read in the existing velocity distributions and create the interpolation functions that will be called when computing DM rate

class VelDist:
    def __init__(self, dist, rho):
        """
        Requires a dist name defined in the velocity distributions files
        Assume that we have both h and g distributions for a given distribution
        """
        self.dist = dist
        self.rho = rho # DM density [GeV/cm3]
        self.gfile= np.loadtxt(f'./velocity_distributions/g{self.dist}.dat')
        self.hfile= np.loadtxt(f'./velocity_distributions/h{self.dist}.dat')

    def gdist(self,vm):
        """
        Interpolate given data to find the integral value for the min velocity vm
        g integral files are typically given in units of [s/km] but we require our values to be unitless so include appropriate c factors
        """
        return np.interp(vm, self.gfile[:, 0], self.gfile[:, 1], right=0) * c * 1E-3
    
    def hdist(self,vm):
        """
        Interpolate given data to find the integral value for the min velocity vm
        h integral files are typically given in units of [km/s] but we require our values to be unitless so include appropriate c factors
        """
        return np.interp(vm, self.hfile[:, 0], self.hfile[:, 1], right=0) / (c * 1E-3)

    def vesc_km_s(self):
        try:
            return vesc_km_s[self.dist]
        except Exception:
            raise RuntimeError(f'Could not find escape velocity for {self.dist}')
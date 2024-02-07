import numpy as np
from constants import *

#### Read in the existing velocity distributions and create the interpolation functions that will be called when computing DM rate

class SolarNDist:
    def __init__(self, dist, flux):
        """
        Requires a dist name defined in the velocity distributions files
        Assume that we have both h and g distributions for a given distribution
        """
        self.flux = flux # DM density [GeV/cm3]
        self.file= np.loadtxt("./snvmodels/solarneutrinospectrum/"+dist+".dat")

    def dist(self,E):
        """
        Interpolate given data to find the integral value for the min velocity vm
        g integral files are typically given in units of [s/km] but we require our values to be unitless so include appropriate c factors
        """
        return np.interp(E,self.file[:,0],self.file[:,1],right=0)*self.flux
    
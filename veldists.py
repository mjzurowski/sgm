import numpy as np
c = 299792458. # speed of light [km/s]
#### Read in the existing velocity distributions and create the interpolation functions that will be called when computing DM rate

class VelDist:
    def __init__(self, dist, rho, norm = c*1E-3):
        """
        Requires a dist name defined in the velocity distributions files
        """
        self.rho = rho # DM density [GeV/cm3]
        self.file= np.loadtxt("./velocity_distributions/"+dist+".dat")
        self.norm=norm # integral files are typically given in units of [s/km] but we require our values to be unitless

    def dist(self,vm):
        """
        Interpolate given data to find the integral value for the min velocity vm
        """
        return np.interp(vm,self.file[:,0],self.file[:,1],right=0)*self.norm
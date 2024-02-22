import numpy as np
from detector import Detector
from targets.na import Na
from targets.i import I
from constants import *

class SABRE(Detector):
    def Nuclei(self):
        return [Na(), I()]
    
    def ER_E(self, E):
        return [E * keV / 0.3, E * keV / 0.09]
    
    def dERdE(self, E):
        return [1. / 0.3, 1. / 0.09]
    
    def Emax(self):
        return 20.
    
    def DeltaE(self, E):
        return np.sqrt(E / 1000.) * 0.014 * 1000.
    
    def Res(self, E1, E2):
        A = 1. / (np.sqrt(2. * np.pi) * self.DeltaE(E2))
        return A * np.exp(-0.5 * pow((E1 - E2) / self.DeltaE(E2), 2.))
    
    def Eff(self, E):
        return np.ones_like(E)
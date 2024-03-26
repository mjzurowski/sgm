import numpy as np
from LAB_detector import LABDetector
from constants import *

import os

class SABREVeto(LABDetector):
    def Emax(self):
        return 20.
    
    def DeltaE(self, E):
        return np.sqrt(E) * 0.154
    
    def Res(self, E1, E2):
        A = 1. / (np.sqrt(2. * np.pi) * self.DeltaE(E2))
        return A * np.exp(-0.5 * pow((E1 - E2) / self.DeltaE(E2), 2.))
    
    def Eff(self, E, coincidence='2fold'):
        if coincidence == '2fold':
            path = os.path.join(os.path.dirname(__file__), 'effs_2fold.npz')
        elif coincidence == '1fold':
            path = os.path.join(os.path.dirname(__file__), 'effs_1fold.npz')
        effs_data = np.load(path)

        energies = effs_data['energies']
        effs = effs_data['effs'] / 100.

        eff_interp = np.interp(E, energies, effs)
        eff_interp = np.where(E > 0.25, 1., eff_interp)

        return eff_interp
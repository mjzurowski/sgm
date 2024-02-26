import numpy as np
from constants import *

class PseudoscalarBSW():    
    def dRdER(self, Target, ER, mX, g_ae):
        C_a = 1.2e19
        dsigdER = C_a / Target.A() * g_ae**2 * (mX  / keV) * (Target.sigma_PE(ER)) * cm2_to_barn
        
        return dsigdER
    
class VectorBSW():    
    def dRdER(self, Target, ER, mX, kappa):
        C_V = 4e23
        dsigdER = C_V / Target.A() * kappa**2 * (keV / mX) * (Target.sigma_PE(ER)) * cm2_to_barn
        
        return dsigdER

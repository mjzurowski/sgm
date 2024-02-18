from dmmodel import DMModel
import numpy as np 
from constants import *

class SIWIMPBSW2(DMModel):
    def dRdER(self, Target, ER, mX, sig, VelDist):
        if ER!=mX:
            return 0
        else:
            C_a = 1.2 * 10**(19)
            g_ae = sig
            dsigdER = (C_a/(Target.A))*(g_ae**2)*mX*(Target.sigma_PE(ER))
            return dsigdER
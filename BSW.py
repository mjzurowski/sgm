from dmmodel import DMModel
import numpy as np 
import scipy.signal 
from target import Target
from targets.na import Na
from targets.i import I 
from constants import *
from BSWDet import DAMA_BSW

class SIWIMPBSW(DMModel):
    def create_delta_function(self,length, places):
        delta_function= np.zeros(length)
        for idx in places:
            delta_func += scipy.signal.unit_impulse(length,idx)
        return delta_function
    def dRdER_a(self, Target, ER, mX, sig, VelDist):
        C_a = 1.2 * 10**(19)
        g_ae = sig
        peaks = [200,700]
        dsigdER = (C_a/(Na.A+I.A))*(g_ae**2)*mX*(Na.sigma_PE+I.sigma_PE)*self.create_delta_function(DAMA_BSW.ROI,peaks)
        return dsigdER
    def dRdER_v(self, Target, ER, mX, sig, VelDist):
        C_v = 4 * (10**23)
        alpha_v = sig
        peaks_v = [200,700]
        dsigdER = (C_v/Na.A+I.A)*(137*alpha_v) *(1/mX)*(Na.sigma_PE+I.sigma_PE)*self.create_delta_function(DAMA_BSW.ROI,peaks_v)
        return dsigdER
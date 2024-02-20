from dmmodel import DMModel
import numpy as np 
import scipy.signal 
from targets.na import Na
from targets.i import I 
from constants import *

class SIWIMPBSW(DMModel):
    def ROI(self):
        return [10,1000]
    def create_delta_function(self,length, places):
        delta_function= np.zeros(length)
        for idx in places:
            delta_func += scipy.signal.unit_impulse(length,idx)
        return delta_function
    def dRdER_a(self,mX, sig):
        C_a = 1.2 * 10**(19)
        g_ae = sig
        peaks = [200,700]
        dsigdER = (C_a/(Na.A+I.A))*(g_ae**2)*mX*(Na.sigma_PE+I.sigma_PE)*self.create_delta_function(self.ROI,peaks)
        return dsigdER
    def dRdER_v(self,mX,sig):
        C_v = 4 * (10**23)
        alpha_v = sig
        peaks_v = [200,700]
        dsigdER = (C_v/Na.A+I.A)*(137*alpha_v)*(1/mX)*(Na.sigma_PE+I.sigma_PE)*self.create_delta_function(self.ROI,peaks_v)
        return dsigdER
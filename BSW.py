from dmmodel import DMModel
import numpy as np 
import scipy.signal 
from constants import *

class SIWIMPBSW(DMModel):
    def ROI(self):
        return [10,1000]
    def create_delta_function(self,length, places):
        delta_function= np.zeros(length)
        for idx in places:
            delta_func += scipy.signal.unit_impulse(length,idx)
        return delta_function
    def dRdER_a(self,Target,ER,mX,sig,VelDist):
        C_a = 1.2 * 10**(19)
        g_ae = sig
        peaks_a = [200,700]
        dsigdER = (C_a/(Target.A))*(g_ae**2)*mX*(Target.sigma_PE(ER))*self.create_delta_function(self.ROI,peaks_a)
        return dsigdER
    def dRdER_v(self,Target,ER,mX,sig,VelDist):
        C_v = 4 * (10**23)
        alpha_v = sig
        peaks_v = [200,700]
        dsigdER = (C_v/Target.A)*(137*alpha_v)*(1/mX)*(Target.sigma_PE(ER))*self.create_delta_function(self.ROI,peaks_v)
        return dsigdER
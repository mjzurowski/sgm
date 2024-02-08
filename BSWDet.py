import numpy as np
from detector import Detector
from targets.na import Na
from targets.i import I
from constants import *
import scipy.integrate as integrate


class DAMA_BSW(Detector):
    def Nuclei(self):
        return[Na(),I()]
    
    def ER_E(self,E):
        #Quenching Factor

        return [E*keV/0.3,E*keV/0.09]
    
    def dERdE(self, E):
        return [1/0.3,1/0.09]
    
    def dRdE_True(self,E,Func,NR=True,**kwargs):
        """
        Rate as a function of observed energy
        Assume "Func" object is something like a background or signal model that is a function of E and target. (NB background will probably be independent of target so it will be a useless variable there)
        Output units are the same as Func
        """ 
        if len(self.ER_E(E))!=len(self.Nuclei()):
            print("You haven't identified an energy transformation for all the target nuclei.")
            ## add some kind of system exit function.
        TotalRate = 0
        if NR:
            for i in range(0,len(self.Nuclei())):
                ER = self.ER_E(E)[i] # energy conversion for this nucleus
                T = self.Nuclei()[i] # target object for computing DM rate
                dERdE = self.dERdE(E)[i]
                TotalRate+=dERdE*float(Func(T,ER,**kwargs))*T.mT()/self.DetMass()
        else:
            for i in range(0,len(self.Nuclei())):
                T = self.Nuclei()[i] # target object for computing DM rate
                TotalRate+=float(Func(T,ER,**kwargs))*T.mT()/self.DetMass()
        return TotalRate
    
    def ROI(self):
        return [10,1000]
    
    def Emax(self):
        return 700
    
    def DeltaE(self,E):
        return (0.488*pow(E,0.5))+(0.0091*E)
    
    def Res(self,E1,E2):
        A = 1./(np.sqrt(2.*np.pi)*self.DeltaE(E2))
        R = A*np.exp(-0.5*pow((E1 - E2)/self.DeltaE(E2), 2.))
        return R
    def dRdE(self,E,Model,NR=True,**kwargs):
        integral = integrate.quad(lambda E2 :self.dRdE_True(E2,Model,NR,**kwargs)*self.Res(E,E2),0,2*self.Emax(),points=self.ROI(),limit=int(1E8))[0]
        return integral 
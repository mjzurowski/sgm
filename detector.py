from abc import ABC, abstractmethod
import scipy as sp
import scipy.integrate as integrate

class Detector(ABC):
### Target/Nucleus specific terms
    ### Need to think of a better way to zip a bunch of target nuclei objects together to then loop over in rate.
    @abstractmethod
    def Nuclei(self):
        """
        Give a list of the target nuclei that make up the detector
        This should be a list of Target objects.
        """
        pass

    def DetMass(self):
        """
        Total detector mass (based on nuclei chosen) [molar mass]
        """
        mass = 0
        for N in self.Nuclei():
            mass+=N.mT()
        return mass

    @abstractmethod
    def ER_E(self,E):
        """
        Recoil energy as a function of observed energy
        (e.g., quenching factor, ionisation yield)
        This should be a list the same length as Targets 
        Output units are [eV] 
        Input units are [keV_0]
        """
        pass

    @abstractmethod
    def dERdE(self,E):
        """
        Derivative of recoil energy wrt observed energy
        Output units should be [keV]/[keV_0] where keV_0 are the observed energy units
        """
        pass

    def dRdE_True(self,Func,E):
        """
        Rate as a function of observed energy
        Assume "Func" object is something like a background or signal model that is a function of E and target. (NB background will probably be independent of target so it will be a useless variable there)
        Output units are the same as Func
        """ 
        if len(self.ER_E(E))!=len(self.Nuclei()):
            print("You haven't identified an energy transformation for all the target nuclei.")
            ## add some kind of system exit function.
        TotalRate = 0
        for i in range(0,len(self.Nuclei())):
            ER = self.ER_E(E)[i] # energy conversion for this nucleus
            T = self.Nuclei()[i] # target object for computing DM rate
            dERdE = self.dERdE(E)[i]
            TotalRate+=dERdE*float(Func(T,ER))*T.mT()/self.DetMass()
        return TotalRate

### General detector terms
    @abstractmethod
    def ROI(self):
        """
        Energy region of interest ([keV_0]) to try and help the integral focus
        """
        pass

    @abstractmethod
    def Emax(self):
        """
        Max energy of operation (used to define range of integrals)
        """
        pass

    @abstractmethod
    def Eff(self,E):
        """
        Efficiency as a function of observed energy
        """
        pass

    @abstractmethod
    def DeltaE(self,E):
        """
        Resolution as a function of observed energy
        """
        pass

    @abstractmethod
    def Res(self,E1,E2):
        """
        Resolution smearing function (probably a gaussian but lets leave arbitrary for now...)
        """
        pass

    def dRdE(self,Model,E):
        """
        Observed rate for Model object smeared with resolution
        """
        integral = integrate.quad(lambda E2: self.dRdE_True(Model,E2)*self.Res(E,E2),0,2*self.Emax(),points=self.ROI(),limit=int(1E8))[0] ## this integral could probs be optimised
        return integral*self.Eff(E)


    
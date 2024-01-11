from abc import ABC, abstractmethod

c = 299792458. # speed of light in m/s
cpd_conversion = 86400*c*1E14 # for conversion to units of cpd/kg/keV

class DMModel(ABC):
    @abstractmethod
    def vmin(self,Target,mX,ER):
        """
        Minimum velocity that can produce recoil energy ER in [km/s]
        Target: target nucleus
        mX: DM mass [eV]
        ER: recoil energy [eV]
        """
        pass

    @abstractmethod
    def dsigdER(self,Target,mX,ER,sig):
        """
        Differential cross section as a function of recoil energy [cm^2]/[eV] 
        Inputs:
            Target: target nucleus (see target.py)
            mX: DM mass [eV]
            ER: recoil energy [eV]
            sig: DM cross section [cm]^2
        """
        pass
    
    def dRdER(self,Target,mX,ER,sig,VelDist):
        """
        Interaction rate as a function of recoil energy in counts/[day]/[kg]/[keV]
        Inputs:
            Target: target nucleus (see target.py)
            mX: DM mass [eV]
            ER: recoil energy [eV]
            sig: DM cross section [cm]^2
            dist: velocity distribution [unitless]
            rho: DM density [GeV]/[cm]^3
        """
        vm = self.vmin(Target,mX,ER)
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*self.dsigdER(Target,mX,ER,sig)*VelDist.dist(vm) ### Need to add normalisation terms

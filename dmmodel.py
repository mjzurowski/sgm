from abc import ABC, abstractmethod
from constants import *

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
        return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*self.dsigdER(Target,mX,ER,sig)*VelDist.dist(vm)
    ## need to do something where we can use the g and h forms...
        # if len(self.dsigdER(Target,mX,ER,sig))!=len(VelDist.dist(vm)):
        #     print("Incorrect definition/combination of velocity distribution and form factors.")
        #     return 0
        # if len(self.dsigdER(Target,mX,ER,sig))==1:
        #     # dependence on only g velocity integral
        #     return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*self.dsigdER(Target,mX,ER,sig)*VelDist.dist(vm)
        # else:
        #     # need to sum g and h contributions
        #     #### I don't love this implementation, but can't think of anything better right now...
        #     return cpd_conversion*Target.N_T()*(VelDist.rho/mX)*(self.dsigdER(Target,mX,ER,sig)[0]*VelDist.dist(vm)[0]+self.dsigdER(Target,mX,ER,sig)[1]*VelDist.dist(vm)[1])


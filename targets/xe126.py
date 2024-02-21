from target import Target
import numpy as np

##### Load useful files (e.g., electron transition probabilities)


#### Define target functions

class Xe126(Target):
    def A(self):
        return 126

    def Z(self):
        return 54

    def FMpp(self,ER):
        pass

    def FMnp(self,ER):
        pass

    def FMpn(self,ER):
        pass

    def FMnn(self,ER):
        pass

    def FS1pp(self,ER):
        pass

    def FS1np(self,ER):
        pass

    def FS1pn(self,ER):
        pass

    def FS1nn(self,ER):
        pass

    def FS2pp(self,ER):
        pass

    def FS2np(self,ER):
        pass

    def FS2pn(self,ER):
        pass

    def FS2nn(self,ER):
        pass

    def FDpp(self,ER):
        pass

    def FDnp(self,ER):
        pass

    def FDpn(self,ER):
        pass

    def FDnn(self,ER):
        pass

    def FS1Dpp(self,ER):
        pass

    def FS1Dnp(self,ER):
        pass

    def FS1Dpn(self,ER):
        pass

    def FS1Dnn(self,ER):
        pass

    def eNL(self):
        """
        List of transition energies for the target
        Should have units of eV
        """
        pass

    def eTransE_E(self):
        """
        Files that give electron transition probability as a function of the kinetic energy of the emitted electron E_e (E_e = E_EM - nuclear bonding)
        Note that targets will have multiple of these depending on the transition energies allowed
        These should have units of [eV]^-1 and be passed as a list the same length as eNL
        """
        pass

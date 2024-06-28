from target import Target
import numpy as np

##### Load useful files (e.g., electron transition probabilities)

# n_10 = np.loadtxt("./targets/elec_prob/na_nl10.dat")
# n_20 = np.loadtxt("./targets/elec_prob/na_nl20.dat")
# n_21 = np.loadtxt("./targets/elec_prob/na_nl21.dat")
# n_30 = np.loadtxt("./targets/elec_prob/na_nl30.dat")

# na_pe = np.loadtxt("./targets/pe_abs/na_pe.dat")

#### Define target functions
## Silicon is made up of a number of different isotopes: 28Si (92.2%), 29Si (4.7%), 30Si (3.1%).
## Lets define them all separately, then create just "Si" to sum them all together in the right ratio
## Credit to Raghda Abdel Khaleq for the detailed shell models.

class Si28(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 28

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1pp(self, ER):
        return 0
    
    def FS1pn(self, ER):
        return 0

    def FS1np(self, ER):
        return 0

    def FS1nn(self, ER):
        return 0

    def FS2pp(self, ER):
        return 0
    
    def FS2pn(self, ER):
        return 0

    def FS2np(self, ER):
        return 0

    def FS2nn(self, ER):
        return 0

    def FPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FPhipp(self, ER):
        return 0
    
    def FPhipn(self, ER):
        return 0

    def FPhinp(self, ER):
        return 0

    def FPhinn(self, ER):
        return 0

    def FDpp(self, ER):
        return 0
    
    def FDpn(self, ER):
        return 0

    def FDnp(self, ER):
        return 0

    def FDnn(self, ER):
        return 0

    def FMPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1Dpp(self, ER):
        return 0

    def FS1Dpn(self, ER):
        return 0

    def FS1Dnp(self, ER):
        return 0

    def FS1Dnn(self, ER):
        return 0   

    ## Will define these three later...
    def eNL(self):
        pass

    def eTransE_E(self):
        pass
    
    def sigma_PE(self,E_gam):
        pass


####################
class Si29(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 29

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(195.994 -335.987*y+195.357*pow(y,2.)-44.0261 *pow(y,3.)+3.36526 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(195.998 -335.997*y+195.406 *pow(y,2.)-44.0637 *pow(y,3.)+3.37088 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(209.992 -366.649*y+219.759 *pow(y,2.)-52.1025 *pow(y,3.)+4.22607 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(209.996 -366.658*y+219.796 *pow(y,2.)-52.1321 *pow(y,3.)+4.23088 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(209.992 -366.649*y+219.759 *pow(y,2.)-52.1025 *pow(y,3.)+4.22607 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(209.996 -366.658*y+219.796 *pow(y,2.)-52.1321 *pow(y,3.)+4.23088 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(224.99 -399.978*y+246.876 *pow(y,2.)-61.4301 *pow(y,3.)+5.30706 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(224.992 -399.984*y+246.9 *pow(y,2.)-61.4494 *pow(y,3.)+5.3103 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)* (0.0000331056 -0.000163055*y+0.000199284 *pow(y,2.)+3.66554e-6 *pow(y,3.)+1.67306e-8 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.00204833 -0.00444303*y+0.00359312 *pow(y,2.)-0.00128387 *pow(y,3.)+0.000171036 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FS1pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.0021637+0.0096352*y-0.0133306 *pow(y,2.)+0.00673258 *pow(y,3.)+0.0000623427 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.0198569 -0.0588155*y+0.0694612 *pow(y,2.)-0.0360334 *pow(y,3.)+0.00673044 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.0021637+0.0096352*y-0.0133306 *pow(y,2.)+0.00673258 *pow(y,3.)+0.0000623427 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.0198569 -0.0588155*y+0.0694612 *pow(y,2.)-0.0360334 *pow(y,3.)+0.00673044 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.141414 -0.56296*y+0.922776 *pow(y,2.)-0.721542 *pow(y,3.)+0.232306 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.192497 -0.722794*y+1.13008 *pow(y,2.)-0.847818 *pow(y,3.)+0.26485 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0000165528 +0.0000306321*y-0.0000519747 *pow(y,2.)-0.0000612043 *pow(y,3.)+0.0000660817 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.00102417 -0.00375031*y+0.00485202 *pow(y,2.)-0.00259768 *pow(y,3.)+0.000491371 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FS2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.00108185-0.000980401*y+0.00134241 *pow(y,2.)-0.000816813 *pow(y,3.)+0.00167485 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.00992846 -0.0206121*y+0.0187405 *pow(y,2.)-0.0152473 *pow(y,3.)+0.00513046 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.00108185-0.000980401*y+0.00134241 *pow(y,2.)-0.000816813 *pow(y,3.)+0.00167485 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.00992846 -0.0206121*y+0.0187405 *pow(y,2.)-0.0152473 *pow(y,3.)+0.00513046 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0707069 -0.0026951*y+0.109597 *pow(y,2.)-0.00208824 *pow(y,3.)+0.0424495 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.0962483 -0.0471923*y+0.149393 *pow(y,2.)-0.0352067 *pow(y,3.)+0.0535677 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(7.52516 -6.02013*y+1.20403 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(8.0029 -6.40235*y+1.28048 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(7.62376 -6.09901*y+1.2198 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(8.17622 -6.541*y+1.30821 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(7.62376 -6.09901*y+1.2198 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(8.17622 -6.541*y+1.30821 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(7.72367 -6.17894*y+1.23579 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(8.35329 -6.68266*y+1.33654 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FPhipp(self, ER):
        return 0
    
    def FPhipn(self, ER):
        return 0

    def FPhinp(self, ER):
        return 0

    def FPhinn(self, ER):
        return 0

    def FDpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.000176249 -0.000140999*y+0.0000281999 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(1.26457e-6 -1.01166e-6*y+2.02332e-7 * pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FDpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.00328811 -0.00263049*y+0.000526098 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.00026013 -0.000208104*y+0.0000416209 * pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FDnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.00328811 -0.00263049*y+0.000526098 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.00026013 -0.000208104* y+0.0000416209* pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FDnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0613433 -0.0490746*y+0.00981492 *pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(0.0535104 -0.0428083 * y+0.00856167* pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-38.4042+48.2793*y-18.1994 *pow(y,2.)+2.01292 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-39.6061+49.7905*y-18.7733 *pow(y,2.)+2.07764 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-38.9075+48.912*y-18.4378 *pow(y,2.)+2.0393 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-40.4632+50.868*y-19.1796 *pow(y,2.)+2.1226 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-41.1471+53.0337*y-20.9495 *pow(y,2.)+2.52781 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-42.4351+54.6941*y-21.6074 *pow(y,2.)+2.60776 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-41.6863+53.7286*y-21.224 *pow(y,2.)+2.56094 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-43.3534+55.8777*y-22.075 *pow(y,2.)+2.6642 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1Dpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0000763861 -0.000218666*y+0.0000735275 *pow(y,2.)+6.86877e-7 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-0.0000504142+0.000074823*y-0.000036418 *pow(y,2.)+5.82204e-6 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1Dpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.00142506 -0.00407945*y+0.00137173 *pow(y,2.)+0.0000128144 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-0.0104582+0.0155218*y-0.00755479 *pow(y,2.)+0.00120776 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1Dnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.0049924+0.0119342*y-0.0103736 *pow(y,2.)+0.00255949 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-0.000489559+0.00111492*y-0.0009418 *pow(y,2.)+0.000229665 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1Dnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.0931385+0.222645*y-0.193531 *pow(y,2.)+0.04775 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-0.101557+0.231286*y-0.195373 *pow(y,2.)+0.0476431 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0  

    ## Will define these three later...
    def eNL(self):
        pass

    def eTransE_E(self):
        pass
    
    def sigma_PE(self,E_gam):
        pass


####################
class Si30(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 30

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(195.996 -335.991*y+194.701 *pow(y,2.)-43.4623 *pow(y,3.)+3.27957 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(196.001 -336.007*y+194.784 *pow(y,2.)-43.5249 *pow(y,3.)+3.2888 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(223.993 -397.318*y+242.474 *pow(y,2.)-58.709 *pow(y,3.)+4.8518 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(223.995 -397.326*y+242.532 *pow(y,2.)-58.7566 *pow(y,3.)+4.86 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(223.993 -397.318*y+242.474 *pow(y,2.)-58.709 *pow(y,3.)+4.8518 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(223.995 -397.326*y+242.532 *pow(y,2.)-58.7566 *pow(y,3.)+4.86 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(255.99 -469.311*y+300.829 *pow(y,2.)-78.5857 *pow(y,3.)+7.17776 *pow(y,4.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(255.986 -469.307*y+300.853 *pow(y,2.)-78.608 *pow(y,3.)+7.18185 *pow(y,4.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1pp(self, ER):
        return 0
    
    def FS1pn(self, ER):
        return 0

    def FS1np(self, ER):
        return 0

    def FS1nn(self, ER):
        return 0

    def FS2pp(self, ER):
        return 0
    
    def FS2pn(self, ER):
        return 0

    def FS2np(self, ER):
        return 0

    def FS2nn(self, ER):
        return 0

    def FPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(8.61265 -6.89012* y+1.37802* pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(9.13079 -7.30459* y+1.46091 *pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(7.07826-5.66261* y+1.13252* pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(7.58577 -6.06858* y+1.21371 * pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(7.07826 -5.66261 * y+1.13252 * pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(7.58577 -6.06858 *y+1.21371 * pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(5.81724 -4.65379 *y+0.930757* pow(y,2.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(6.30219 -5.04171 * y+1.00834 * pow(y,2.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
        
    def FPhipp(self, ER):
        return 0
    
    def FPhipn(self, ER):
        return 0

    def FPhinp(self, ER):
        return 0

    def FPhinn(self, ER):
        return 0

    def FDpp(self, ER):
        return 0
    
    def FDpn(self, ER):
        return 0

    def FDnp(self, ER):
        return 0

    def FDnn(self, ER):
        return 0

    def FMPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-41.0858+51.6505*y-19.4011 *pow(y,2.)+2.12587 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-42.3041+53.1822*y-19.9837 *pow(y,2.)+2.19178 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-33.7662+42.4487*y-15.9447 *pow(y,2.)+1.74713 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-35.1452+44.1825*y-16.602 *pow(y,2.)+1.82088 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-46.9548+61.8234*y-25.0791 *pow(y,2.)+3.14501 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-48.3476+63.6577*y-25.8258 *pow(y,2.)+3.23935 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-38.5896+50.8093*y-20.6112 *pow(y,2.)+2.58471 *pow(y,3.))
        elif self.shell_model=="USDB":
            return np.exp(-2.*y)*(-40.166+52.8853*y-21.4555 *pow(y,2.)+2.69117 *pow(y,3.))
        else:
            print(self.shell_model+" doesn't exist for Si.")
            return 0

    def FS1Dpp(self, ER):
        return 0

    def FS1Dpn(self, ER):
        return 0

    def FS1Dnp(self, ER):
        return 0

    def FS1Dnn(self, ER):
        return 0   

    ## Will define these three later...
    def eNL(self):
        pass

    def eTransE_E(self):
        pass
    
    def sigma_PE(self,E_gam):
        pass


####################
class Si(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        # Initialise all the isotopes: 70Ge (20%), 72Ge (27.4%), 73Ge (7.76%), 74Ge (36.5%), 76Ge (7.75%).
        # We'll assume we always want to use the same shell model, but somewhat trivial to tweak if we decide to change that
        self.si28 = Si28(shell_model)
        self.si29 = Si29(shell_model)
        self.si30 = Si30(shell_model)

    def A(self):
        return 0.922*self.si28.A()+0.047*self.si29.A()+0.031*self.si30.A()

    def FMpp(self, ER):
        return 0.922*self.si28.FMpp(ER)+0.047*self.si29.FMpp(ER)+0.031*self.si30.FMpp(ER)

    def FMpn(self, ER):
        return 0.922*self.si28.FMpn(ER)+0.047*self.si29.FMpn(ER)+0.031*self.si30.FMpn(ER)
        
    def FMnp(self, ER):
        return 0.922*self.si28.FMnp(ER)+0.047*self.si29.FMnp(ER)+0.031*self.si30.FMnp(ER)
        
    def FMnn(self, ER):
        return 0.922*self.si28.FMnn(ER)+0.047*self.si29.FMnn(ER)+0.031*self.si30.FMnn(ER)

    def FS1pp(self, ER):
        return 0.922*self.si28.FS1pp(ER)+0.047*self.si29.FS1pp(ER)+0.031*self.si30.FS1pp(ER)
    
    def FS1pn(self, ER):
        return 0.922*self.si28.FS1pn(ER)+0.047*self.si29.FS1pn(ER)+0.031*self.si30.FS1pn(ER)

    def FS1np(self, ER):
        return 0.922*self.si28.FS1np(ER)+0.047*self.si29.FS1np(ER)+0.031*self.si30.FS1np(ER)

    def FS1nn(self, ER):
        return 0.922*self.si28.FS1nn(ER)+0.047*self.si29.FS1nn(ER)+0.031*self.si30.FS1nn(ER)

    def FS2pp(self, ER):
        return 0.922*self.si28.FS2pp(ER)+0.047*self.si29.FS2pp(ER)+0.031*self.si30.FS2pp(ER)
    
    def FS2pn(self, ER):
        return 0.922*self.si28.FS2pn(ER)+0.047*self.si29.FS2pn(ER)+0.031*self.si30.FS2pn(ER)

    def FS2np(self, ER):
        return 0.922*self.si28.FS2np(ER)+0.047*self.si29.FS2np(ER)+0.031*self.si30.FS2np(ER)

    def FS2nn(self, ER):
        return 0.922*self.si28.FS2nn(ER)+0.047*self.si29.FS2nn(ER)+0.031*self.si30.FS2nn(ER)

    def FPhi2pp(self, ER):
        return 0.922*self.si28.FPhi2pp(ER)+0.047*self.si29.FPhi2pp(ER)+0.031*self.si30.FPhi2pp(ER)
    
    def FPhi2pn(self, ER):
        return 0.922*self.si28.FPhi2pn(ER)+0.047*self.si29.FPhi2pn(ER)+0.031*self.si30.FPhi2pn(ER)

    def FPhi2np(self, ER):
        return 0.922*self.si28.FPhi2np(ER)+0.047*self.si29.FPhi2np(ER)+0.031*self.si30.FPhi2np(ER)

    def FPhi2nn(self, ER):
        return 0.922*self.si28.FPhi2nn(ER)+0.047*self.si29.FPhi2nn(ER)+0.031*self.si30.FPhi2nn(ER)
        
    def FPhipp(self, ER):
        return 0.922*self.si28.FPhipp(ER)+0.047*self.si29.FPhipp(ER)+0.031*self.si30.FPhipp(ER)
    
    def FPhipn(self, ER):
        return 0.922*self.si28.FPhipn(ER)+0.047*self.si29.FPhipn(ER)+0.031*self.si30.FPhipn(ER)

    def FPhinp(self, ER):
        return 0.922*self.si28.FPhinp(ER)+0.047*self.si29.FPhinp(ER)+0.031*self.si30.FPhinp(ER)

    def FPhinn(self, ER):
        return 0.922*self.si28.FPhinn(ER)+0.047*self.si29.FPhinn(ER)+0.031*self.si30.FPhinn(ER)

    def FDpp(self, ER):
        return 0.922*self.si28.FDpp(ER)+0.047*self.si29.FDpp(ER)+0.031*self.si30.FDpp(ER)
    
    def FDpn(self, ER):
        return 0.922*self.si28.FDpn(ER)+0.047*self.si29.FDpn(ER)+0.031*self.si30.FDpn(ER)

    def FDnp(self, ER):
        return 0.922*self.si28.FDnp(ER)+0.047*self.si29.FDnp(ER)+0.031*self.si30.FDnp(ER)

    def FDnn(self, ER):
        return 0.922*self.si28.FDnn(ER)+0.047*self.si29.FDnn(ER)+0.031*self.si30.FDnn(ER)

    def FMPhi2pp(self, ER):
        return 0.922*self.si28.FMPhi2pp(ER)+0.047*self.si29.FMPhi2pp(ER)+0.031*self.si30.FMPhi2pp(ER)
    
    def FMPhi2pn(self, ER):
        return 0.922*self.si28.FMPhi2pn(ER)+0.047*self.si29.FMPhi2pn(ER)+0.031*self.si30.FMPhi2pn(ER)

    def FMPhi2np(self, ER):
        return 0.922*self.si28.FMPhi2np(ER)+0.047*self.si29.FMPhi2np(ER)+0.031*self.si30.FMPhi2np(ER)

    def FMPhi2nn(self, ER):
        return 0.922*self.si28.FMPhi2nn(ER)+0.047*self.si29.FMPhi2nn(ER)+0.031*self.si30.FMPhi2nn(ER)

    def FS1Dpp(self, ER):
        return 0.922*self.si28.FS1Dpp(ER)+0.047*self.si29.FS1Dpp(ER)+0.031*self.si30.FS1Dpp(ER)

    def FS1Dpn(self, ER):
        return 0.922*self.si28.FS1Dpn(ER)+0.047*self.si29.FS1Dpn(ER)+0.031*self.si30.FS1Dpn(ER)

    def FS1Dnp(self, ER):
        return 0.922*self.si28.FS1Dnp(ER)+0.047*self.si29.FS1Dnp(ER)+0.031*self.si30.FS1Dnp(ER)

    def FS1Dnn(self, ER):
        return 0.922*self.si28.FS1Dnn(ER)+0.047*self.si29.FS1Dnn(ER)+0.031*self.si30.FS1Dnn(ER) 

## Will define these three later...
    def eNL(self):
        pass

    def eTransE_E(self):
        ## This is going to be potentially annoying to define with multi-isotopes and the way its being done in the back end... need to think about how to deal with ith
        pass
    
    def sigma_PE(self,E_gam):
        pass
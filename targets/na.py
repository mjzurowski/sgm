from target import Target
import numpy as np

##### Load useful files (e.g., electron transition probabilities)

n_10 = np.loadtxt("./targets/elec_prob/na_nl10.dat")
n_20 = np.loadtxt("./targets/elec_prob/na_nl20.dat")
n_21 = np.loadtxt("./targets/elec_prob/na_nl21.dat")
n_30 = np.loadtxt("./targets/elec_prob/na_nl30.dat")

#### Define target functions

class Na(Target):
    def A(self):
        return 23

    def Z(self):
        return 11

    def FMpp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(120.-180.*y+87.*np.power(y,2.)-17.*np.power(y,3.)+1.2*np.power(y,4.))

    def FMnp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(130.-200.*y+100.*np.power(y,2.)-20.*np.power(y,3.)+1.5*np.power(y,4.))

    def FMpn(self,ER):
        return self.FMnp(ER)

    def FMnn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(140.-220.*y+120.*np.power(y,2.)-25.*np.power(y,3.)+1.8*np.power(y,4.))

    def FS1pp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.273-0.824*y+1.19*np.power(y,2.)-0.477*np.power(y,3.)+0.0593*np.power(y,4.))

    def FS1np(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.0219-0.0578*y+0.036*np.power(y,2.)-0.003*np.power(y,3.)+0.000363*np.power(y,4.))

    def FS1pn(self,ER):
        return self.FS1np(ER)

    def FS1nn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.00176-0.00396*y+0.00228*np.power(y,2.)-0.0000195*np.power(y,3.))

    def FS2pp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.136-0.267*y+0.458*np.power(y,2.)-0.112*np.power(y,3.)+0.00828*np.power(y,4.))

    def FS2np(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.011-0.03*y+0.0217*np.power(y,2.)-0.00897*np.power(y,3.)+0.000592*np.power(y,4.))

    def FS2pn(self,ER):
        return self.FS2np(ER)

    def FS2nn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.000882-0.0031*y+0.00399*np.power(y,2.)-0.00203*np.power(y,3.)+0.000409*np.power(y,4.))

    def FDpp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.231-0.185*y+0.0502*np.power(y,2.))

    def FDnp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.0812-0.065*y+0.0138*np.power(y,2.))

    def FDpn(self,ER):
        return self.FDnp(ER)

    def FDnn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.0286-0.0228*y+0.00462*np.power(y,2.))

    def FS1Dpp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.25+0.48*y-0.29*np.power(y,2.)+0.049*np.power(y,3.))

    def FS1Dnp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.02+0.031*y-0.0076*np.power(y,2.)+0.00027*np.power(y,3.))

    def FS1Dpn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.088+0.17*y-0.081*np.power(y,2.)+0.011*np.power(y,3.))

    def FS1Dnn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.0071+0.011*y-0.003*np.power(y,2.))
    
    def eNL(self):
        return [1.1E3,6.6E1,3.9E1,6.1]

    def eTransE_E(self):
        return [n_10,n_20,n_21,n_30]
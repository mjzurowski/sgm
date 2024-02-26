from target import Target
import numpy as np

##### Load useful files (e.g., electron transition probabilities)
i_10 = np.loadtxt("./targets/elec_prob/i_nl10.dat")
i_20 = np.loadtxt("./targets/elec_prob/i_nl20.dat")
i_21 = np.loadtxt("./targets/elec_prob/i_nl21.dat")
i_30 = np.loadtxt("./targets/elec_prob/i_nl30.dat")
i_31 = np.loadtxt("./targets/elec_prob/i_nl31.dat")
i_32 = np.loadtxt("./targets/elec_prob/i_nl32.dat")
i_40 = np.loadtxt("./targets/elec_prob/i_nl40.dat")
i_41 = np.loadtxt("./targets/elec_prob/i_nl41.dat")
i_42 = np.loadtxt("./targets/elec_prob/i_nl42.dat")
i_50 = np.loadtxt("./targets/elec_prob/i_nl50.dat")
i_51 = np.loadtxt("./targets/elec_prob/i_nl51.dat")

i_pe = np.loadtxt("./targets/pe_abs/i_pe.dat")

##### Define target functions

class I(Target):
    def A(self):
        return 127

    def FMpp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(2800.-10000.*y+14000.*pow(y,2.)-9800.*pow(y,3.)+3800.*pow(y,4.)-840.*pow(y,5.)+100.*pow(y,6.)-6.3*pow(y,7.)+0.15*pow(y,8.))

    def FMnp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(3900.-15000.*y+23000.*pow(y,2.)-18000.*pow(y,3.)+7900.*pow(y,4.)-2000.*pow(y,5.)+290.*pow(y,6.)-23.*pow(y,7.)+0.75*pow(y,8.)-0.0048*pow(y,9.))

    def FMpn(self,ER):
        return self.FMnp(ER)

    def FMnn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(5500.-23000.*y+38000.*pow(y,2.)-32000.*pow(y,3.)+16000.*pow(y,4.)-4600.*pow(y,5.)+790.*pow(y,6.)-75.*pow(y,7.)+3.3*pow(y,8.)-0.041*pow(y,9.)+0.00015*pow(y,10.))

    def FS1pp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.26-1.6*y+5.3*pow(y,2.)-8.9*pow(y,3.)+8.7*pow(y,4.)-4.9*pow(y,5.)+1.5*pow(y,6.)-0.25*pow(y,7.)+0.016*pow(y,8.))

    def FS1np(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.065-0.46*y+1.3*pow(y,2.)-1.8*pow(y,3.)+1.4*pow(y,4.)-0.65*pow(y,5.)+0.17*pow(y,6.)-0.026*pow(y,7.)+0.002*pow(y,8.))

    def FS1pn(self,ER):
        return self.FS1np(ER)

    def FS1nn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.016-0.13*y+0.37*pow(y,2.)-0.48*pow(y,3.)+0.34*pow(y,4.)-0.14*pow(y,5.)+0.033*pow(y,6.)-0.0048*pow(y,7.)+0.00041*pow(y,8.))

    def FS2pp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.13-0.49*y+1.8*pow(y,2.)-2.8*pow(y,3.)+2.7*pow(y,4.)-1.6*pow(y,5.)+0.53*pow(y,6.)-0.092*pow(y,7.)+0.0067*pow(y,8.))

    def FS2np(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.032-0.13*y+0.26*pow(y,2.)-0.3*pow(y,3.)+0.21*pow(y,4.)-0.098*pow(y,5.)+0.027*pow(y,6.)-0.0042*pow(y,7.)+0.00033*pow(y,8.))

    def FS2pn(self,ER):
        return self.FS2np(ER)

    def FS2nn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.008-0.032*y+0.053*pow(y,2.)-0.046*pow(y,3.)+0.025*pow(y,4.)-0.0086*pow(y,5.)+0.0019*pow(y,6.)-0.00026*pow(y,7.))

    def FDpp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.54-1.3*y+1.6*pow(y,2.)-1.2*pow(y,3.)+0.51*pow(y,4.)-0.11*pow(y,5.)+0.0097*pow(y,6.))

    def FDnp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.23-0.65*y+0.79*pow(y,2.)-0.54*pow(y,3.)+0.2*pow(y,4.)-0.04*pow(y,5.)+0.0039*pow(y,6.)-0.00014*pow(y,7.))

    def FDpn(self,ER):
        return self.FDnp(ER)

    def FDnn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(0.1-0.32*y+0.4*pow(y,2.)-0.25*pow(y,3.)+0.084*pow(y,4.)-0.016*pow(y,5.)+0.0018*pow(y,6.)-0.00011*pow(y,7.))

    def FS1Dpp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.37+1.6*y-3.1*pow(y,2.)+3.4*pow(y,3.)-2.1*pow(y,4.)+0.73*pow(y,5.)-0.13*pow(y,6.)+0.0086*pow(y,7.))

    def FS1Dnp(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.093+0.48*y-0.85*pow(y,2.)+0.79*pow(y,3.)-0.43*pow(y,4.)+0.13*pow(y,5.)-0.021*pow(y,6.)+0.0017*pow(y,7.))

    def FS1Dpn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.16+0.75*y-1.5*pow(y,2.)+1.5*pow(y,3.)-0.82*pow(y,4.)+0.26*pow(y,5.)-0.047*pow(y,6.)+0.0043*pow(y,7.)-0.00015*pow(y,8.))

    def FS1Dnn(self,ER):
        y = self.Y(ER)
        return np.exp(-2.*y)*(-0.04+0.22*y-0.43*pow(y,2.)+0.38*pow(y,3.)-0.19*pow(y,4.)+0.053*pow(y,5.)-0.009*pow(y,6.)+0.00088*pow(y,7.))
    
    def eNL(self):
        return [3.3E4,5.1E3,4.6E3,1.0E3,8.7E2,6.1E2,1.8E2,1.3E2,5.1E1,1.9E1,8.8]

    def eTransE_E(self):
        return [i_10,i_20,i_21,i_30,i_31,i_32,i_40,i_41,i_42,i_50,i_51]
    
    def sigma_PE(self,E_gam):
        return np.interp(E_gam,i_pe[:,0],i_pe[:,1])

    
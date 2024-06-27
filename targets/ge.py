from target import Target
import numpy as np

##### Load useful files (e.g., electron transition probabilities)

# n_10 = np.loadtxt("./targets/elec_prob/na_nl10.dat")
# n_20 = np.loadtxt("./targets/elec_prob/na_nl20.dat")
# n_21 = np.loadtxt("./targets/elec_prob/na_nl21.dat")
# n_30 = np.loadtxt("./targets/elec_prob/na_nl30.dat")

# na_pe = np.loadtxt("./targets/pe_abs/na_pe.dat")

#### Define target functions
## Germanium is made up of a number of different isotopes: 70Ge (20%), 72Ge (27.4%), 73Ge (7.76%), 74Ge (36.5%), 76Ge (7.75%).
## Lets define them all separately, then create just "Ge" to sum them all together in the right ratio
## Credit to Raghda Abdel Khaleq for the detailed shell models.

class Ge70(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 70

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1023.86 -2818. *y+2913.01 *pow(y,2.)-1428.53 *pow(y,3.)+353.01 *pow(y,4.)-42.0142 *pow(y,5.)+1.92717 *pow(y,6.)-0.00270658 *pow(y,7.)+9.65284e-7 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1023.41 -2825.68 *y+2951.99 *pow(y,2.)-1478.77 *pow(y,3.)+378.034 *pow(y,4.)-47.4449 *pow(y,5.)+2.40216 *pow(y,6.)-0.0139098 *pow(y,7.)+0.0000214293 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1023.29 -2828.58 *y+2926.66 *pow(y,2.)-1431.06 *pow(y,3.)+352.403 *pow(y,4.)-42.1794 *pow(y,5.)+2.05998 *pow(y,6.)-0.0163005 *pow(y,7.)+0.0000353456 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1215.49 -3492.33 *y+3796.45 *pow(y,2.)-1980.47 *pow(y,3.)+528.605 *pow(y,4.)-69.751 *pow(y,5.)+3.80204 *pow(y,6.)-0.0337146 *pow(y,7.)+0.0000223007 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1214.79 -3502.95 *y+3837.97 *pow(y,2.)-2030.76 *pow(y,3.)+554.24 *pow(y,4.)-75.7766 *pow(y,5.)+4.42145 *pow(y,6.)-0.0550292 *pow(y,7.)+0.000132836 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1214.42 -3523.21 *y+3865.2 *pow(y,2.)-2036.23 *pow(y,3.)+551.766 *pow(y,4.)-75.2609 *pow(y,5.)+4.53602 *pow(y,6.)-0.0767866 *pow(y,7.)+0.000262912 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1215.49 -3492.33 *y+3796.45 *pow(y,2.)-1980.47 *pow(y,3.)+528.605 *pow(y,4.)-69.751 *pow(y,5.)+3.80204 *pow(y,6.)-0.0337146 *pow(y,7.)+0.0000223007 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1214.79 -3502.95 *y+3837.97 *pow(y,2.)-2030.76 *pow(y,3.)+554.24 *pow(y,4.)-75.7766 *pow(y,5.)+4.42145 *pow(y,6.)-0.0550292 *pow(y,7.)+0.000132836 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1214.42 -3523.21 *y+3865.2 *pow(y,2.)-2036.23 *pow(y,3.)+551.766 *pow(y,4.)-75.2609 *pow(y,5.)+4.53602 *pow(y,6.)-0.0767866 *pow(y,7.)+0.000262912 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1442.98 -4320.37 *y+4929.61 *pow(y,2.)-2728.02 *pow(y,3.)+783.524 *pow(y,4.)-113.896 *pow(y,5.)+7.2311 *pow(y,6.)-0.113198 *pow(y,7.)+0.000515206 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1441.95 -4334.7 *y+4973.73 *pow(y,2.)-2774.86 *pow(y,3.)+806.609 *pow(y,4.)-119.613 *pow(y,5.)+7.92387 *pow(y,6.)-0.147742 *pow(y,7.)+0.000823428 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1441.24 -4378.64 *y+5079.28 *pow(y,2.)-2870.21 *pow(y,3.)+850.322 *pow(y,4.)-130.675 *pow(y,5.)+9.43336 *pow(y,6.)-0.240444 *pow(y,7.)+0.00195562 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(35.7749 -57.61 *y+31.6045 *pow(y,2.)-6.79615 *pow(y,3.)+0.513349 *pow(y,4.)-0.00276336 *pow(y,5.)+3.86114e-6 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(56.5726 -92.7096 *y+54.0115 *pow(y,2.)-13.2732 *pow(y,3.)+1.24951 *pow(y,4.)-0.0197305 *pow(y,5.)+0.0000857171 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(42.4481 -70.3572 *y+41.2081 *pow(y,2.)-10.1447 *pow(y,3.)+0.984154 *pow(y,4.)-0.0219989 *pow(y,5.)+0.000141383 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(25.3832 -45.0209 *y+28.5792 *pow(y,2.)-7.84106 *pow(y,3.)+0.908725 *pow(y,4.)-0.0338264 *pow(y,5.)+0.0000892027 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(39.9244 -71.4517 *y+46.6257 *pow(y,2.)-13.5302 *pow(y,3.)+1.70871 *pow(y,4.)-0.0725119 *pow(y,5.)+0.000531345 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(47.8393 -86.9937 *y+57.4149 *pow(y,2.)-16.8203 *pow(y,3.)+2.18029 *pow(y,4.)-0.102593 *pow(y,5.)+0.00105165 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(25.3832 -45.0209 *y+28.5792 *pow(y,2.)-7.84106 *pow(y,3.)+0.908725 *pow(y,4.)-0.0338264 *pow(y,5.)+0.0000892027 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(39.9244 -71.4517 *y+46.6257 *pow(y,2.)-13.5302 *pow(y,3.)+1.70871 *pow(y,4.)-0.0725119 *pow(y,5.)+0.000531345 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(47.8393 -86.9937 *y+57.4149 *pow(y,2.)-16.8203 *pow(y,3.)+2.18029 *pow(y,4.)-0.102593 *pow(y,5.)+0.00105165 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(18.01 -34.8846 *y+25.1251 *pow(y,2.)-8.3584 *pow(y,3.)+1.31397 *pow(y,4.)-0.0880645 *pow(y,5.)+0.00206082 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(28.1754 -54.6766 *y+39.551 *pow(y,2.)-13.2472 *pow(y,3.)+2.09646 *pow(y,4.)-0.140826 *pow(y,5.)+0.00329371 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(53.9153 -106.721 *y+78.4707 *pow(y,2.)-26.6939 *pow(y,3.)+4.33836 *pow(y,4.)-0.30907 *pow(y,5.)+0.00782248 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(-191.386+417.477 *y-325.597 *pow(y,2.)+112.561 *pow(y,3.)-17.4277 *pow(y,4.)+1.00317 *pow(y,5.)-0.00339742 *pow(y,6.)+1.93057e-6 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-240.618+529.337 *y-424.008 *pow(y,2.)+155.128 *pow(y,3.)-26.3827 *pow(y,4.)+1.77435 *pow(y,5.)-0.0188424 *pow(y,6.)+0.0000428585 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-208.415+460.773 *y-367.293 *pow(y,2.)+132.241 *pow(y,3.)-22.0203 *pow(y,4.)+1.48091 *pow(y,5.)-0.0218002 *pow(y,6.)+0.0000706913 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-135.793+318.386 *y-276.609 *pow(y,2.)+112.563 *pow(y,3.)-22.4268 *pow(y,4.)+2.03098 *pow(y,5.)-0.0634823 *pow(y,6.)+0.0000446014 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-169.809+399.188 *y-349.799 *pow(y,2.)+144.617 *pow(y,3.)-29.5025 *pow(y,4.)+2.76552 *pow(y,5.)-0.091904 *pow(y,6.)+0.000265672 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-234.885+557.104 *y-488.74 *pow(y,2.)+200.549 *pow(y,3.)-40.461 *pow(y,4.)+3.78216 *pow(y,5.)-0.131636 *pow(y,6.)+0.000525824 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-227.206+523.074 *y-434.08 *pow(y,2.)+162.468 *pow(y,3.)-27.9509 *pow(y,4.)+1.90652 *pow(y,5.)-0.02086 *pow(y,6.)+0.0000446014 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-285.613+663.324 *y-562.175 *pow(y,2.)+219.789 *pow(y,3.)-40.6866 *pow(y,4.)+3.1291 *pow(y,5.)-0.0544103 *pow(y,6.)+0.000265672 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-247.342+580.708 *y-496.972 *pow(y,2.)+196.215 *pow(y,3.)-37.0177 *pow(y,4.)+3.02826 *pow(y,5.)-0.0732338 *pow(y,6.)+0.000525824 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-161.208+397.46 *y-365.295 *pow(y,2.)+159.203 *pow(y,3.)-34.5762 *pow(y,4.)+3.52518 *pow(y,5.)-0.135214 *pow(y,6.)+0.00103041 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-201.563+498.537 *y-460.49 *pow(y,2.)+202.246 *pow(y,3.)-44.4093 *pow(y,4.)+4.60301 *pow(y,5.)-0.182948 *pow(y,6.)+0.00164686 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-278.756+699.333 *y-655.005 *pow(y,2.)+291.921 *pow(y,3.)-65.5354 *pow(y,4.)+7.11409 *pow(y,5.)-0.317711 *pow(y,6.)+0.00391124 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
class Ge72(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 72

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1023.98 -2821.39 *y+2965.33 *pow(y,2.)-1509.59 *pow(y,3.)+395.325 *pow(y,4.)-50.985 *pow(y,5.)+2.59908 *pow(y,6.)-0.00687813 *pow(y,7.)+4.6747e-6 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1023.39 -2825.57 *y+2938.34 *pow(y,2.)-1456.2 *pow(y,3.)+366.126 *pow(y,4.)-44.9447 *pow(y,5.)+2.22207 *pow(y,6.)-0.0133092 *pow(y,7.)+0.000021295 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1023.31 -2826.64 *y+2917.67 *pow(y,2.)-1419.67 *pow(y,3.)+346.819 *pow(y,4.)-40.9926 *pow(y,5.)+1.95851 *pow(y,6.)-0.0138169 *pow(y,7.)+0.0000264679 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1279.93 -3692.09 *y+4070.13 *pow(y,2.)-2179.85 *pow(y,3.)+602.572 *pow(y,4.)-82.408 *pow(y,5.)+4.49675 *pow(y,6.)-0.0171605 *pow(y,7.)+0.0000153234 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1278.64 -3744.31 *y+4160.43 *pow(y,2.)-2230.4 *pow(y,3.)+617.476 *pow(y,4.)-86.2208 *pow(y,5.)+5.28229 *pow(y,6.)-0.0847285 *pow(y,7.)+0.000225984 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1278.17 -3762.37 *y+4187.08 *pow(y,2.)-2238.58 *pow(y,3.)+616.741 *pow(y,4.)-85.986 *pow(y,5.)+5.38365 *pow(y,6.)-0.102263 *pow(y,7.)+0.00032888 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1279.93 -3692.09 *y+4070.13 *pow(y,2.)-2179.85 *pow(y,3.)+602.572 *pow(y,4.)-82.408 *pow(y,5.)+4.49675 *pow(y,6.)-0.0171605 *pow(y,7.)+0.0000153234 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1278.64 -3744.31 *y+4160.43 *pow(y,2.)-2230.4 *pow(y,3.)+617.476 *pow(y,4.)-86.2208 *pow(y,5.)+5.28229 *pow(y,6.)-0.0847285 *pow(y,7.)+0.000225984 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1278.17 -3762.37 *y+4187.08 *pow(y,2.)-2238.58 *pow(y,3.)+616.741 *pow(y,4.)-85.986 *pow(y,5.)+5.38365 *pow(y,6.)-0.102263 *pow(y,7.)+0.00032888 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1599.85 -4821.79 *y+5568.71 *pow(y,2.)-3134.7 *pow(y,3.)+914.288 *pow(y,4.)-132.629 *pow(y,5.)+7.75785 *pow(y,6.)-0.0385976 *pow(y,7.)+0.0000502291 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1597.55 -4945.57 *y+5854.1 *pow(y,2.)-3381.26 *pow(y,3.)+1024.93 *pow(y,4.)-161.083 *pow(y,5.)+11.8312 *pow(y,6.)-0.299456 *pow(y,7.)+0.00239816 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1596.49 -4988.87 *y+5960.3 *pow(y,2.)-3478.19 *pow(y,3.)+1070. *pow(y,4.)-172.765 *pow(y,5.)+13.4874 *pow(y,6.)-0.408067 *pow(y,7.)+0.00408655 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(67.762 -109.541 *y+63.6637 *pow(y,2.)-15.7471 *pow(y,3.)+1.44526 *pow(y,4.)-0.010188 *pow(y,5.)+0.0000186988 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(50.7782 -83.3168 *y+48.8247 *pow(y,2.)-12.1489 *pow(y,3.)+1.16432 *pow(y,4.)-0.0189721 *pow(y,5.)+0.0000851801 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(39.3855 -65.0509 *y+37.9757 *pow(y,2.)-9.30851 *pow(y,3.)+0.890904 *pow(y,4.)-0.0182241 *pow(y,5.)+0.000105871 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(6.92791 -12.9797 *y+8.79139 *pow(y,2.)-2.65951 *pow(y,3.)+0.360765 *pow(y,4.)-0.0176616 *pow(y,5.)+0.0000612935 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(43.8843 -82.1023 *y+56.0174 *pow(y,2.)-17.1401 *pow(y,3.)+2.34401 *pow(y,4.)-0.116037 *pow(y,5.)+0.000903937 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(51.2546 -95.9683 *y+65.2955 *pow(y,2.)-19.8231 *pow(y,3.)+2.69294 *pow(y,4.)-0.135791 *pow(y,5.)+0.00131552 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(6.92791 -12.9797 *y+8.79139 *pow(y,2.)-2.65951 *pow(y,3.)+0.360765 *pow(y,4.)-0.0176616 *pow(y,5.)+0.0000612935 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(43.8843 -82.1023 *y+56.0174 *pow(y,2.)-17.1401 *pow(y,3.)+2.34401 *pow(y,4.)-0.116037 *pow(y,5.)+0.000903937 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(51.2546 -95.9683 *y+65.2955 *pow(y,2.)-19.8231 *pow(y,3.)+2.69294 *pow(y,4.)-0.135791 *pow(y,5.)+0.00131552 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.708302 -1.50906 *y+1.17896 *pow(y,2.)-0.423531 *pow(y,3.)+0.0750997 *pow(y,4.)-0.00631895 *pow(y,5.)+0.000200916 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(37.9263 -79.6819 *y+62.3648 *pow(y,2.)-22.7544 *pow(y,3.)+4.0408 *pow(y,4.)-0.326226 *pow(y,5.)+0.00959263 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(66.7005 -139.612 *y+108.882 *pow(y,2.)-39.5825 *pow(y,3.)+6.9963 *pow(y,4.)-0.560845 *pow(y,5.)+0.0163462 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(-263.414+575.805 *y-462.45 *pow(y,2.)+171.401 *pow(y,3.)-29.6009 *pow(y,4.)+1.95719 *pow(y,5.)-0.00942512 *pow(y,6.)+9.34941e-6 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-227.96+501.716 *y-401.097 *pow(y,2.)+146.238 *pow(y,3.)-24.7426 *pow(y,4.)+1.65167 *pow(y,5.)-0.0180523 *pow(y,6.)+0.0000425901 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)* (-200.758+443.061 *y-352.034 *pow(y,2.)+126.111 *pow(y,3.)-20.8135 *pow(y,4.)+1.37116 *pow(y,5.)-0.018373 *pow(y,6.)+0.0000529357 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-26.9312+65.7907 *y-60.094 *pow(y,2.)+25.9336 *pow(y,3.)-5.61175 *pow(y,4.)+0.582807 *pow(y,5.)-0.023028 *pow(y,6.)+0.0000306468 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-197.011+478.93 *y-434.079 *pow(y,2.)+185.463 *pow(y,3.)-39.4004 *pow(y,4.)+3.9439 *pow(y,5.)-0.148924 *pow(y,6.)+0.000451968 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-261.257+634.25 *y-571.067 *pow(y,2.)+240.975 *pow(y,3.)-50.2749 *pow(y,4.)+4.91908 *pow(y,5.)-0.182968 *pow(y,6.)+0.000657761 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-329.255+762.3 *y-647.339 *pow(y,2.)+254.584 *pow(y,3.)-46.9403 *pow(y,4.)+3.35955 *pow(y,5.)-0.0201238 *pow(y,6.)+0.0000306468 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-284.817+674.521 *y-583.412 *pow(y,2.)+233.951 *pow(y,3.)-44.851 *pow(y,4.)+3.6628 *pow(y,5.)-0.0785517 *pow(y,6.)+0.000451968 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-250.756+598.873 *y-520.941 *pow(y,2.)+209.515 *pow(y,3.)-40.4457 *pow(y,4.)+3.42343 *pow(y,5.)-0.0894524 *pow(y,6.)+0.000657761 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-33.6627+86.5877 *y-83.3181 *pow(y,2.)+37.9867 *pow(y,3.)-8.69494 *pow(y,4.)+0.956283 *pow(y,5.)-0.0401773 *pow(y,6.)+0.000100458 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-246.149+639.579 *y-622.929 *pow(y,2.)+289.785 *pow(y,3.)-68.3614 *pow(y,4.)+7.89174 *pow(y,5.)-0.381013 *pow(y,6.)+0.00479631 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-326.323+851.38 *y-832.066 *pow(y,2.)+388.746 *pow(y,3.)-92.4036 *pow(y,4.)+10.8473 *pow(y,5.)-0.548278 *pow(y,6.)+0.0081731 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
class Ge73(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 73

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1023.96 -2820.72 *y+2943.43 *pow(y,2.)-1474.31 *pow(y,3.)+376.997 *pow(y,4.)-47.2283 *pow(y,5.)+2.34033 *pow(y,6.)-0.00581515 *pow(y,7.)+3.73787e-6 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1022.57 -2821.85 *y+2926.09 *pow(y,2.)-1441.23 *pow(y,3.)+358.989 *pow(y,4.)-43.5153 *pow(y,5.)+2.1176 *pow(y,6.)-0.0119992 *pow(y,7.)+0.000018122 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1022.02 -2820.04 *y+2905.24 *pow(y,2.)-1408.94 *pow(y,3.)+342.429 *pow(y,4.)-40.1322 *pow(y,5.)+1.88306 *pow(y,6.)-0.0112592 *pow(y,7.)+0.0000181689 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1311.8 -3823.75 *y+4247.43 *pow(y,2.)-2285.42 *pow(y,3.)+634.832 *pow(y,4.)-88.0003 *pow(y,5.)+5.08999 *pow(y,6.)-0.0495837 *pow(y,7.)+0.00005318 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1307.9 -3855.48 *y+4309.43 *pow(y,2.)-2322.47 *pow(y,3.)+646.467 *pow(y,4.)-90.9714 *pow(y,5.)+5.67141 *pow(y,6.)-0.0980863 *pow(y,7.)+0.00025333 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1306.4 -3870.91 *y+4334.98 *pow(y,2.)-2332.17 *pow(y,3.)+646.979 *pow(y,4.)-91.032 *pow(y,5.)+5.79116 *pow(y,6.)-0.114791 *pow(y,7.)+0.000321909 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1311.8 -3823.75 *y+4247.43 *pow(y,2.)-2285.42 *pow(y,3.)+634.832 *pow(y,4.)-88.0003 *pow(y,5.)+5.08999 *pow(y,6.)-0.0495837 *pow(y,7.)+0.00005318 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1307.9 -3855.48 *y+4309.43 *pow(y,2.)-2322.47 *pow(y,3.)+646.467 *pow(y,4.)-90.9714 *pow(y,5.)+5.67141 *pow(y,6.)-0.0980863 *pow(y,7.)+0.00025333 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1306.4 -3870.91 *y+4334.98 *pow(y,2.)-2332.17 *pow(y,3.)+646.979 *pow(y,4.)-91.032 *pow(y,5.)+5.79116 *pow(y,6.)-0.114791 *pow(y,7.)+0.000321909 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1680.56 -5167.82 *y+6095.35 *pow(y,2.)-3513.6 *pow(y,3.)+1057.64 *pow(y,4.)-161.794 *pow(y,5.)+10.8041 *pow(y,6.)-0.174453 *pow(y,7.)+0.00122931 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1672.85 -5246.25 *y+6296.25 *pow(y,2.)-3692.64 *pow(y,3.)+1139.94 *pow(y,4.)-183.643 *pow(y,5.)+14.0588 *pow(y,6.)-0.393753 *pow(y,7.)+0.00361486 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1669.91 -5288.25 *y+6404.82 *pow(y,2.)-3792.59 *pow(y,3.)+1186.84 *pow(y,4.)-196.06 *pow(y,5.)+15.8817 *pow(y,6.)-0.519048 *pow(y,7.)+0.00572743 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS1pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.000204083 -0.000463558 *y+0.00148348 *pow(y,2.)-0.00229056 *pow(y,3.)+0.00205049 *pow(y,4.)-0.000756444 *pow(y,5.)+0.000102527 *pow(y,6.)-2.05037e-6 *pow(y,7.)+3.00769e-8 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.000178729 -0.000879632 *y+0.00262286 *pow(y,2.)-0.00411546 *pow(y,3.)+0.00414324 *pow(y,4.)-0.00165035 *pow(y,5.)+0.000258082 *pow(y,6.)-0.0000104406 *pow(y,7.)+1.30346e-7 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.0000175837 -0.0000904842 *y+0.000558419 *pow(y,2.)-0.00125076 *pow(y,3.)+0.00305584 *pow(y,4.)-0.00160437 *pow(y,5.)+0.000270483 *pow(y,6.)-9.82312e-6 *pow(y,7.)+1.05683e-7 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FS1pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0122577 -0.0531362 *y+0.0723673 *pow(y,2.)-0.0669153 *pow(y,3.)+0.038403 *pow(y,4.)-0.0119103 *pow(y,5.)+0.00184147 *pow(y,6.)-0.000139208 *pow(y,7.)+4.97489e-6 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.011069 -0.0626324 *y+0.163422 *pow(y,2.)-0.232366 *pow(y,3.)+0.156508 *pow(y,4.)-0.0508548 *pow(y,5.)+0.00778239 *pow(y,6.)-0.000497255 *pow(y,7.)+0.0000107912 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.00318903 -0.0183393 *y+0.0742487 *pow(y,2.)-0.157254 *pow(y,3.)+0.123229 *pow(y,4.)-0.04216 *pow(y,5.)+0.00653537 *pow(y,6.)-0.000416954 *pow(y,7.)+8.04581e-6 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS1np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0122577 -0.0531362 *y+0.0723673 *pow(y,2.)-0.0669153 *pow(y,3.)+0.038403 *pow(y,4.)-0.0119103 *pow(y,5.)+0.00184147 *pow(y,6.)-0.000139208 *pow(y,7.)+4.97489e-6 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.011069 -0.0626324 *y+0.163422 *pow(y,2.)-0.232366 *pow(y,3.)+0.156508 *pow(y,4.)-0.0508548 *pow(y,5.)+0.00778239 *pow(y,6.)-0.000497255 *pow(y,7.)+0.0000107912 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.00318903 -0.0183393 *y+0.0742487 *pow(y,2.)-0.157254 *pow(y,3.)+0.123229 *pow(y,4.)-0.04216 *pow(y,5.)+0.00653537 *pow(y,6.)-0.000416954 *pow(y,7.)+8.04581e-6 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS1nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.736223 -4.71068 *y+11.6726 *pow(y,2.)-12.5199 *pow(y,3.)+6.91299 *pow(y,4.)-2.0712 *pow(y,5.)+0.347982 *pow(y,6.)-0.030887 *pow(y,7.)+0.00191149 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.68552 -4.38401 *y+10.9223 *pow(y,2.)-11.797 *pow(y,3.)+6.60331 *pow(y,4.)-2.0122 *pow(y,5.)+0.348628 *pow(y,6.)-0.032337 *pow(y,7.)+0.00199453 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.578373 -3.6759 *y+9.02308 *pow(y,2.)-9.6299 *pow(y,3.)+5.35495 *pow(y,4.)-1.63303 *pow(y,5.)+0.284934 *pow(y,6.)-0.0267136 *pow(y,7.)+0.00167647 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.000102041 -0.00096119 *y+0.00335189 *pow(y,2.)-0.00364089 *pow(y,3.)+0.0019991 *pow(y,4.)-0.000462685 *pow(y,5.)+0.000040244 *pow(y,6.)-3.45152e-7 *pow(y,7.)+1.2863e-9 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.0000893644 -0.000566113 *y+0.0013045 *pow(y,2.)-0.00143037 *pow(y,3.)+0.000884122 *pow(y,4.)-0.000259919 *pow(y,5.)+0.0000309428 *pow(y,6.)-7.51951e-7 *pow(y,7.)+7.34639e-9 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(8.79184e-6-0.000128723 *y+0.000563505 *pow(y,2.)-0.000701525 *pow(y,3.)+0.000434789 *pow(y,4.)-0.000115875 *pow(y,5.)+0.000013001 *pow(y,6.)-3.4661e-7 *pow(y,7.)+3.46682e-9 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FS2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.00612884 -0.0386733 *y+0.0575045 *pow(y,2.)-0.0535467 *pow(y,3.)+0.0257882 *pow(y,4.)-0.0065118 *pow(y,5.)+0.000816404 *pow(y,6.)-0.0000402374 *pow(y,7.)+7.36917e-7 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.00553449 -0.0264412 *y+0.0489084 *pow(y,2.)-0.0505379 *pow(y,3.)+0.0266135 *pow(y,4.)-0.00720555 *pow(y,5.)+0.000953058 *pow(y,6.)-0.0000519915 *pow(y,7.)+1.91813e-6 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.00159452 -0.0144716 *y+0.0313208 *pow(y,2.)-0.0339063 *pow(y,3.)+0.0173419 *pow(y,4.)-0.00457872 *pow(y,5.)+0.000606076 *pow(y,6.)-0.0000333043 *pow(y,7.)+1.00314e-6 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.00612884 -0.0386733 *y+0.0575045 *pow(y,2.)-0.0535467 *pow(y,3.)+0.0257882 *pow(y,4.)-0.0065118 *pow(y,5.)+0.000816404 *pow(y,6.)-0.0000402374 *pow(y,7.)+7.36917e-7 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.00553449 -0.0264412 *y+0.0489084 *pow(y,2.)-0.0505379 *pow(y,3.)+0.0266135 *pow(y,4.)-0.00720555 *pow(y,5.)+0.000953058 *pow(y,6.)-0.0000519915 *pow(y,7.)+1.91813e-6 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.00159452 -0.0144716 *y+0.0313208 *pow(y,2.)-0.0339063 *pow(y,3.)+0.0173419 *pow(y,4.)-0.00457872 *pow(y,5.)+0.000606076 *pow(y,6.)-0.0000333043 *pow(y,7.)+1.00314e-6 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.368112 -1.17814 *y+2.27237 *pow(y,2.)-1.98163 *pow(y,3.)+1.01909 *pow(y,4.)-0.296483 *pow(y,5.)+0.056608 *pow(y,6.)-0.00599025 *pow(y,7.)+0.000951084 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.34276 -1.10375 *y+2.16199 *pow(y,2.)-1.90745 *pow(y,3.)+1.03078 *pow(y,4.)-0.314771 *pow(y,5.)+0.0649442 *pow(y,6.)-0.00731164 *pow(y,7.)+0.000989408 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.289187 -1.0152 *y+1.96819 *pow(y,2.)-1.7061 *pow(y,3.)+0.884432 *pow(y,4.)-0.260285 *pow(y,5.)+0.0522085 *pow(y,6.)-0.00580833 *pow(y,7.)+0.000824121 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(51.127 -82.7519 *y+47.2943 *pow(y,2.)-11.1964 *pow(y,3.)+0.969881 *pow(y,4.)-0.00726876 *pow(y,5.)+0.0000149515 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(46.5381 -76.2937 *y+44.6246 *pow(y,2.)-11.0456 *pow(y,3.)+1.04839 *pow(y,4.)-0.0164757 *pow(y,5.)+0.0000724881 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(36.9372 -60.7344 *y+35.1576 *pow(y,2.)-8.48263 *pow(y,3.)+0.788995 *pow(y,4.)-0.0142748 *pow(y,5.)+0.0000726757 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(17.696 -34.1357 *y+23.579 *pow(y,2.)-7.21471 *pow(y,3.)+0.987163 *pow(y,4.)-0.0493724 *pow(y,5.)+0.00021272 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(45.8377 -86.9084 *y+60.018 *pow(y,2.)-18.5882 *pow(y,3.)+2.58122 *pow(y,4.)-0.131055 *pow(y,5.)+0.00101332 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(54.1357 -102.272 *y+70.0974 *pow(y,2.)-21.4037 *pow(y,3.)+2.92081 *pow(y,4.)-0.147507 *pow(y,5.)+0.00128763 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(17.696 -34.1357 *y+23.579 *pow(y,2.)-7.21471 *pow(y,3.)+0.987163 *pow(y,4.)-0.0493724 *pow(y,5.)+0.00021272 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(45.8377 -86.9084 *y+60.018 *pow(y,2.)-18.5882 *pow(y,3.)+2.58122 *pow(y,4.)-0.131055 *pow(y,5.)+0.00101332 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(54.1357 -102.272 *y+70.0974 *pow(y,2.)-21.4037 *pow(y,3.)+2.92081 *pow(y,4.)-0.147507 *pow(y,5.)+0.00128763 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(6.56241 -15.2582 *y+13.3402 *pow(y,2.)-5.47127 *pow(y,3.)+1.13804 *pow(y,4.)-0.115782 *pow(y,5.)+0.00491724 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(45.2289 -97.5041 *y+78.3722 *pow(y,2.)-29.4332 *pow(y,3.)+5.41721 *pow(y,4.)-0.459504 *pow(y,5.)+0.0144594 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(79.3514 -169.382 *y+134.836 *pow(y,2.)-50.1323 *pow(y,3.)+9.10103 *pow(y,4.)-0.754992 *pow(y,5.)+0.0229097 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FPhipp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0189154 -0.0144994 *y+0.00611477 *pow(y,2.)-0.00130759 *pow(y,3.)+0.000135463 *pow(y,4.)+2.4551e-6 *pow(y,5.)+1.51325e-8 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.000502647 -0.00142759 *y+0.00121767 *pow(y,2.)-0.000268569 *pow(y,3.)+0.0000160394 *pow(y,4.)+1.07185e-6 *pow(y,5.)+4.8707e-8 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.000354089 -0.00149899 *y+0.00207648 *pow(y,2.)-0.00085719 *pow(y,3.)+0.000119079 *pow(y,4.)-9.77543e-7 *pow(y,5.)+1.45538e-8 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FPhipn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.10485+0.130524 *y-0.0712295 *pow(y,2.)+0.0213711 *pow(y,3.)-0.00312987 *pow(y,4.)+0.000158349 *pow(y,5.)-2.67567e-7 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-0.0073194+0.0172517 *y-0.0126229 *pow(y,2.)+0.0039925 *pow(y,3.)-0.000571078 *pow(y,4.)+0.0000259874 *pow(y,5.)+3.41694e-7 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-0.00168767+0.0061126 *y-0.00710558 *pow(y,2.)+0.00268132 *pow(y,3.)-0.000387668 *pow(y,4.)+0.0000192225 *pow(y,5.)+1.2328e-7 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhinp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.10485+0.130524 *y-0.0712295 *pow(y,2.)+0.0213711 *pow(y,3.)-0.00312987 *pow(y,4.)+0.000158349 *pow(y,5.)-2.67567e-7 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-0.0073194+0.0172517 *y-0.0126229 *pow(y,2.)+0.0039925 *pow(y,3.)-0.000571078 *pow(y,4.)+0.0000259874 *pow(y,5.)+3.41694e-7 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-0.00168767+0.0061126 *y-0.00710558 *pow(y,2.)+0.00268132 *pow(y,3.)-0.000387668 *pow(y,4.)+0.0000192225 *pow(y,5.)+1.2328e-7 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhinn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.581194 -1.00151 *y+0.867587 *pow(y,2.)-0.361044 *pow(y,3.)+0.0887292 *pow(y,4.)-0.0110701 *pow(y,5.)+0.000860753 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.106583 -0.199717 *y+0.156802 *pow(y,2.)-0.0598136 *pow(y,3.)+0.0158485 *pow(y,4.)-0.00220302 *pow(y,5.)+0.00016324 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.00804383 -0.0242156 *y+0.0239678 *pow(y,2.)-0.00834067 *pow(y,3.)+0.00222349 *pow(y,4.)-0.000336765 *pow(y,5.)+0.0000206655 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FDpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.0068596 -0.0119041 *y+0.00887033 *pow(y,2.)-0.00313383 *pow(y,3.)+0.00047719 *pow(y,4.)-0.0000153398 *pow(y,5.)+1.39822e-7 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.0119721 -0.0216 *y+0.0152189 *pow(y,2.)-0.00485965 *pow(y,3.)+0.000687894 *pow(y,4.)-0.0000353117 *pow(y,5.)+5.92009e-7 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.0243195 -0.0422941 *y+0.0285927 *pow(y,2.)-0.00900814 *pow(y,3.)+0.00122899 *pow(y,4.)-0.000045406 *pow(y,5.)+5.00332e-7 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FDpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.143256 -0.295723 *y+0.248303 *pow(y,2.)-0.0993772 *pow(y,3.)+0.0184168 *pow(y,4.)-0.00131503 *pow(y,5.)+0.0000197052 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.185956 -0.38918 *y+0.321564 *pow(y,2.)-0.12273 *pow(y,3.)+0.0222375 *pow(y,4.)-0.00172173 *pow(y,5.)+0.0000407036 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.260219 -0.534587 *y+0.419885 *pow(y,2.)-0.154654 *pow(y,3.)+0.0268152 *pow(y,4.)-0.00188275 *pow(y,5.)+0.0000332333 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FDnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(0.143256 -0.295723 *y+0.248303 *pow(y,2.)-0.0993772 *pow(y,3.)+0.0184168 *pow(y,4.)-0.00131503 *pow(y,5.)+0.0000197052 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(0.185956 -0.38918 *y+0.321564 *pow(y,2.)-0.12273 *pow(y,3.)+0.0222375 *pow(y,4.)-0.00172173 *pow(y,5.)+0.0000407036 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(0.260219 -0.534587 *y+0.419885 *pow(y,2.)-0.154654 *pow(y,3.)+0.0268152 *pow(y,4.)-0.00188275 *pow(y,5.)+0.0000332333 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FDnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(2.99178 -7.15989 *y+6.93935 *pow(y,2.)-3.14012 *pow(y,3.)+0.720619 *pow(y,4.)-0.0799055 *pow(y,5.)+0.00357837 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(2.88836 -6.87871 *y+6.6473 *pow(y,2.)-3.00001 *pow(y,3.)+0.692974 *pow(y,4.)-0.0779032 *pow(y,5.)+0.00362051 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(2.78436 -6.59792 *y+6.25673 *pow(y,2.)-2.7726 *pow(y,3.)+0.627064 *pow(y,4.)-0.0690249 *pow(y,5.)+0.00312635 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2pp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-228.777+499.892 *y-395.972 *pow(y,2.)+142.19 *pow(y,3.)-23.3913 *pow(y,4.)+1.45844 *pow(y,5.)-0.00763234 *pow(y,6.)+7.47574e-6 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-218.147+479.775 *y-382.063 *pow(y,2.)+138.218 *pow(y,3.)-23.105 *pow(y,4.)+1.51563 *pow(y,5.)-0.0161181 *pow(y,6.)+0.000036244 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-194.294+427.776 *y-338.299 *pow(y,2.)+120.159 *pow(y,3.)-19.5428 *pow(y,4.)+1.25147 *pow(y,5.)-0.0148279 *pow(y,6.)+0.0000363378 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-79.5167+200.771 *y-189.956 *pow(y,2.)+85.1572 *pow(y,3.)-19.2585 *pow(y,4.)+2.10752 *pow(y,5.)-0.0888734 *pow(y,6.)+0.00010636 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-214.888+528.303 *y-484.058 *pow(y,2.)+208.9 *pow(y,3.)-44.8667 *pow(y,4.)+4.55785 *pow(y,5.)-0.176278 *pow(y,6.)+0.00050666 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-284.769+696.832 *y-632.828 *pow(y,2.)+269.396 *pow(y,3.)-56.7482 *pow(y,4.)+5.61427 *pow(y,5.)-0.210972 *pow(y,6.)+0.000643817 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-293.089+687.417 *y-587.516 *pow(y,2.)+230.794 *pow(y,3.)-42.4386 *pow(y,4.)+3.11741 *pow(y,5.)-0.0349803 *pow(y,6.)+0.00010636 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-279.018+666.177 *y-580.125 *pow(y,2.)+234.168 *pow(y,3.)-45.2637 *pow(y,4.)+3.75612 *pow(y,5.)-0.0854217 *pow(y,6.)+0.00050666 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-248.357+597.413 *y-522.383 *pow(y,2.)+210.92 *pow(y,3.)-40.8612 *pow(y,4.)+3.47971 *pow(y,5.)-0.0923624 *pow(y,6.)+0.000643817 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-101.87+273.206 *y-275.726 *pow(y,2.)+133.047 *pow(y,3.)-32.701 *pow(y,4.)+3.97733 *pow(y,5.)-0.203399 *pow(y,6.)+0.00245862 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-274.85+727.495 *y-722.662 *pow(y,2.)+343.49 *pow(y,3.)-83.1421 *pow(y,4.)+9.94279 *pow(y,5.)-0.50863 *pow(y,6.)+0.00722972 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-364.006+964.883 *y-958.859 *pow(y,2.)+456.185 *pow(y,3.)-110.788 *pow(y,4.)+13.3845 *pow(y,5.)-0.707796 *pow(y,6.)+0.0114549 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS1Dpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.00118319+0.0023704 *y-0.00332808 *pow(y,2.)+0.00230287 *pow(y,3.)-0.000883856 *pow(y,4.)+0.000150282 *pow(y,5.)-5.87253e-6 *pow(y,6.)+6.44573e-8 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-0.00146279+0.00491921 *y-0.00984692 *pow(y,2.)+0.00796406 *pow(y,3.)-0.00279635 *pow(y,4.)+0.000414489 *pow(y,5.)-0.0000195337 *pow(y,6.)+2.74912e-7 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)* (-0.000653931+0.00225116 *y-0.00975989 *pow(y,2.)+0.00980705 *pow(y,3.)-0.00381755 *pow(y,4.)+0.000574597 *pow(y,5.)-0.0000211487 *pow(y,6.)+2.29121e-7 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS1Dpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.0247097+0.0576306 *y-0.0744189 *pow(y,2.)+0.0534968 *pow(y,3.)-0.0202228 *pow(y,4.)+0.00363838 *pow(y,5.)-0.000290423 *pow(y,6.)+9.35121e-6 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-0.0227208+0.0829663 *y-0.1745 *pow(y,2.)+0.162398 *pow(y,3.)-0.0679092 *pow(y,4.)+0.0128812 *pow(y,5.)-0.000973633 *pow(y,6.)+0.0000194848 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-0.00699709+0.0262935 *y-0.112179 *pow(y,2.)+0.135144 *pow(y,3.)-0.0613329 *pow(y,4.)+0.0117098 *pow(y,5.)-0.000855018 *pow(y,6.)+0.0000155952 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS1Dnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-0.0710647+0.289015 *y-0.383118 *pow(y,2.)+0.239052 *pow(y,3.)-0.0746961 *pow(y,4.)+0.0110683 *pow(y,5.)-0.000656255 *pow(y,6.)+9.35121e-6 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-0.090593+0.371402 *y-0.498492 *pow(y,2.)+0.304624 *pow(y,3.)-0.0919462 *pow(y,4.)+0.0136469 *pow(y,5.)-0.000897863 *pow(y,6.)+0.0000194848 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-0.118599+0.480011 *y-0.609646 *pow(y,2.)+0.35663 *pow(y,3.)-0.104884 *pow(y,4.)+0.0150891 *pow(y,5.)-0.000923413 *pow(y,6.)+0.0000155952 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FS1Dnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-1.48412+6.52393 *y-10.0708 *pow(y,2.)+7.11424 *pow(y,3.)-2.58424 *pow(y,4.)+0.494767 *pow(y,5.)-0.0482323 *pow(y,6.)+0.00192315 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)* (-1.40714+6.17499 *y-9.55783 *pow(y,2.)+6.76987 *pow(y,3.)-2.47988 *pow(y,4.)+0.481409 *pow(y,5.)-0.04827 *pow(y,6.)+0.00201417 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-1.26901+5.53621 *y-8.42784 *pow(y,2.)+5.88585 *pow(y,3.)-2.13304 *pow(y,4.)+0.411552 *pow(y,5.)-0.0411436 *pow(y,6.)+0.00171176 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0   

## Will define these three later...
    def eNL(self):
        pass

    def eTransE_E(self):
        pass
    
    def sigma_PE(self,E_gam):
        pass


####################
class Ge74(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 74

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1023.92 -2818.75 *y+2910.8 *pow(y,2.)-1423.58 *pow(y,3.)+350.294 *pow(y,4.)-41.4647 *pow(y,5.)+1.89466 *pow(y,6.)-0.00327916 *pow(y,7.)+1.44671e-6 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1023.66 -2825.66 *y+2919.98 *pow(y,2.)-1425.87 *pow(y,3.)+350.237 *pow(y,4.)-41.629 *pow(y,5.)+1.98107 *pow(y,6.)-0.0116585 *pow(y,7.)+0.0000183561 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1023.63 -2825.75 *y+2911.93 *pow(y,2.)-1412.26 *pow(y,3.)+343.142 *pow(y,4.)-40.1868 *pow(y,5.)+1.88408 *pow(y,6.)-0.0115506 *pow(y,7.)+0.0000190273 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1343.54 -3957.18 *y+4421.02 *pow(y,2.)-2379.94 *pow(y,3.)+659.091 *pow(y,4.)-91.2045 *pow(y,5.)+5.38326 *pow(y,6.)-0.0698114 *pow(y,7.)+0.0000580033 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1343.04 -3995.22 *y+4498.38 *pow(y,2.)-2437.85 *pow(y,3.)+682.1 *pow(y,4.)-96.7638 *pow(y,5.)+6.16924 *pow(y,6.)-0.118462 *pow(y,7.)+0.000324445 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1342.92 -4006.3 *y+4519.47 *pow(y,2.)-2451.36 *pow(y,3.)+686.305 *pow(y,4.)-97.6096 *pow(y,5.)+6.29772 *pow(y,6.)-0.128342 *pow(y,7.)+0.000369663 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1343.54 -3957.18 *y+4421.02 *pow(y,2.)-2379.94 *pow(y,3.)+659.091 *pow(y,4.)-91.2045 *pow(y,5.)+5.38326 *pow(y,6.)-0.0698114 *pow(y,7.)+0.0000580033 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1343.04 -3995.22 *y+4498.38 *pow(y,2.)-2437.85 *pow(y,3.)+682.1 *pow(y,4.)-96.7638 *pow(y,5.)+6.16924 *pow(y,6.)-0.118462 *pow(y,7.)+0.000324445 *pow(y,8.))
        elif self.shell_model=="jj44b":
             return np.exp(-2.*y)*(1342.92 -4006.3 *y+4519.47 *pow(y,2.)-2451.36 *pow(y,3.)+686.305 *pow(y,4.)-97.6096 *pow(y,5.)+6.29772 *pow(y,6.)-0.128342 *pow(y,7.)+0.000369663 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1762.93 -5531.68 *y+6655.74 *pow(y,2.)-3918.75 *pow(y,3.)+1211.35 *pow(y,4.)-193.276 *pow(y,5.)+14.1399 *pow(y,6.)-0.326777 *pow(y,7.)+0.00232555 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1762.07 -5619.54 *y+6858.47 *pow(y,2.)-4094.37 *pow(y,3.)+1290.83 *pow(y,4.)-214.16 *pow(y,5.)+17.2598 *pow(y,6.)-0.545439 *pow(y,7.)+0.0057346 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1761.81 -5648.39 *y+6933.97 *pow(y,2.)-4168.66 *pow(y,3.)+1326.98 *pow(y,4.)-223.564 *pow(y,5.)+18.5498 *pow(y,6.)-0.627129 *pow(y,7.)+0.00718182 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(35.6195 -57.4435 *y+31.8372 *pow(y,2.)-7.0258 *pow(y,3.)+0.551649 *pow(y,4.)-0.00349761 *pow(y,5.)+5.78684e-6 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(41.5175 -68.1671 *y+39.8375 *pow(y,2.)-9.84421 *pow(y,3.)+0.937185 *pow(y,4.)-0.0157678 *pow(y,5.)+0.0000734242 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(36.9847 -60.8467 *y+35.343 *pow(y,2.)-8.59276 *pow(y,3.)+0.806767 *pow(y,4.)-0.0147999 *pow(y,5.)+0.0000761094 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(26.4251 -51.5139 *y+35.8594 *pow(y,2.)-10.9454 *pow(y,3.)+1.48539 *pow(y,4.)-0.0734543 *pow(y,5.)+0.000232013 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(50.8446 -97.7865 *y+68.3577 *pow(y,2.)-21.432 *pow(y,3.)+3.02998 *pow(y,4.)-0.159412 *pow(y,5.)+0.00129778 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(55.9536 -107.024 *y+74.1946 *pow(y,2.)-22.9565 *pow(y,3.)+3.19691 *pow(y,4.)-0.167026 *pow(y,5.)+0.00147865 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(26.4251 -51.5139 *y+35.8594 *pow(y,2.)-10.9454 *pow(y,3.)+1.48539 *pow(y,4.)-0.0734543 *pow(y,5.)+0.000232013 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(50.8446 -97.7865 *y+68.3577 *pow(y,2.)-21.432 *pow(y,3.)+3.02998 *pow(y,4.)-0.159412 *pow(y,5.)+0.00129778 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(55.9536 -107.024 *y+74.1946 *pow(y,2.)-22.9565 *pow(y,3.)+3.19691 *pow(y,4.)-0.167026 *pow(y,5.)+0.00147865 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(19.604 -44.818 *y+37.9067 *pow(y,2.)-14.9041 *pow(y,3.)+2.90288 *pow(y,4.)-0.267743 *pow(y,5.)+0.00930218 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(62.2672 -137.274 *y+112.611 *pow(y,2.)-43.1229 *pow(y,3.)+8.11714 *pow(y,4.)-0.709245 *pow(y,5.)+0.0229384 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(84.6513 -184.564 *y+149.661 *pow(y,2.)-56.6023 *pow(y,3.)+10.5085 *pow(y,4.)-0.903788 *pow(y,5.)+0.0287273 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(-190.975+416.861 *y-325.766 *pow(y,2.)+113.239 *pow(y,3.)-17.7017 *pow(y,4.)+1.03326 *pow(y,5.)-0.00415356 *pow(y,6.)+2.89342e-6 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-206.155+453.771 *y-360.697 *pow(y,2.)+129.857 *pow(y,3.)-21.5506 *pow(y,4.)+1.40439 *pow(y,5.)-0.0156005 *pow(y,6.)+0.0000367121 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-194.573+428.615 *y-339.465 *pow(y,2.)+120.982 *pow(y,3.)-19.7847 *pow(y,4.)+1.27599 *pow(y,5.)-0.0152506 *pow(y,6.)+0.0000380547 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-141.679+356.965 *y-334.501 *pow(y,2.)+147.036 *pow(y,3.)-32.2092 *pow(y,4.)+3.36131 *pow(y,5.)-0.133142 *pow(y,6.)+0.000116007 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-252.469+626.746 *y-578.63 *pow(y,2.)+250.834 *pow(y,3.)-54.0509 *pow(y,4.)+5.51892 *pow(y,5.)-0.216097 *pow(y,6.)+0.00064889 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-294.367+727.203 *y-666.522 *pow(y,2.)+286.103 *pow(y,3.)-60.8796 *pow(y,4.)+6.12135 *pow(y,5.)-0.236035 *pow(y,6.)+0.000739327 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-250.589+595.208 *y-512.169 *pow(y,2.)+200.962 *pow(y,3.)-36.8051 *pow(y,4.)+2.76152 *pow(y,5.)-0.0432081 *pow(y,6.)+0.000116007 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-270.475+653.34 *y-575.206 *pow(y,2.)+234.985 *pow(y,3.)-46.1733 *pow(y,4.)+3.95679 *pow(y,5.)-0.100534 *pow(y,6.)+0.00064889 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-255.264+619.17 *y-546.556 *pow(y,2.)+223.363 *pow(y,3.)-43.9306 *pow(y,4.)+3.81255 *pow(y,5.)-0.104163 *pow(y,6.)+0.000739327 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-185.904+504.168 *y-513.811 *pow(y,2.)+250.097 *pow(y,3.)-62.0033 *pow(y,4.)+7.60737 *pow(y,5.)-0.393713 *pow(y,6.)+0.00465109 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-331.239+893.311 *y-904.024 *pow(y,2.)+437.885 *pow(y,3.)-108.384 *pow(y,4.)+13.3811 *pow(y,5.)-0.72275 *pow(y,6.)+0.0114692 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-386.185+1040.06 *y-1050.55 *pow(y,2.)+508.106 *pow(y,3.)-125.734 *pow(y,4.)+15.5743 *pow(y,5.)-0.853076 *pow(y,6.)+0.0143636 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
class Ge76(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        self.shell_model=shell_model
        
    def A(self):
        return 76

    def FMpp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1023.88 -2817.91 *y+2899.26 *pow(y,2.)-1405.86 *pow(y,3.)+341.22 *pow(y,4.)-39.5982 *pow(y,5.)+1.76111 *pow(y,6.)-0.00242521 *pow(y,7.)+8.48247e-7 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1023.9 -2825.7 *y+2911.06 *pow(y,2.)-1411.34 *pow(y,3.)+342.688 *pow(y,4.)-40.0615 *pow(y,5.)+1.86602 *pow(y,6.)-0.0105226 *pow(y,7.)+0.0000158507 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1023.9 -2824.87 *y+2905.76 *pow(y,2.)-1404.01 *pow(y,3.)+339.038 *pow(y,4.)-39.3042 *pow(y,5.)+1.80704 *pow(y,6.)-0.00953459 *pow(y,7.)+0.0000133881 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMpn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1407.26 -4215.51 *y+4779.75 *pow(y,2.)-2607.9 *pow(y,3.)+733.322 *pow(y,4.)-104.01 *pow(y,5.)+6.52025 *pow(y,6.)-0.112776 *pow(y,7.)+0.0000763214 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1407.65 -4242.36 *y+4832.41 *pow(y,2.)-2646.4 *pow(y,3.)+748.319 *pow(y,4.)-107.621 *pow(y,5.)+7.04276 *pow(y,6.)-0.147115 *pow(y,7.)+0.000393906 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1407.68 -4248.49 *y+4848.53 *pow(y,2.)-2661.45 *pow(y,3.)+754.549 *pow(y,4.)-108.846 *pow(y,5.)+7.15387 *pow(y,6.)-0.151125 *pow(y,7.)+0.000381186 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnp(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1407.26 -4215.51 *y+4779.75 *pow(y,2.)-2607.9 *pow(y,3.)+733.322 *pow(y,4.)-104.01 *pow(y,5.)+6.52025 *pow(y,6.)-0.112776 *pow(y,7.)+0.0000763214 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1407.65 -4242.36 *y+4832.41 *pow(y,2.)-2646.4 *pow(y,3.)+748.319 *pow(y,4.)-107.621 *pow(y,5.)+7.04276 *pow(y,6.)-0.147115 *pow(y,7.)+0.000393906 *pow(y,8.))
        elif self.shell_model=="jj44b":
             return np.exp(-2.*y)*(1407.68 -4248.49 *y+4848.53 *pow(y,2.)-2661.45 *pow(y,3.)+754.549 *pow(y,4.)-108.846 *pow(y,5.)+7.15387 *pow(y,6.)-0.151125 *pow(y,7.)+0.000381186 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
        
    def FMnn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(1934.21 -6264.64 *y+7776.54 *pow(y,2.)-4729.5 *pow(y,3.)+1520.11 *pow(y,4.)-256.886 *pow(y,5.)+20.985 *pow(y,6.)-0.660661 *pow(y,7.)+0.00686705 *pow(y,8.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(1935.23 -6324.02 *y+7909.96 *pow(y,2.)-4844.29 *pow(y,3.)+1571.95 *pow(y,4.)-270.574 *pow(y,5.)+23.0667 *pow(y,6.)-0.813386 *pow(y,7.)+0.009789 *pow(y,8.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(1935.31 -6342.46 *y+7969.43 *pow(y,2.)-4913.96 *pow(y,3.)+1608.89 *pow(y,4.)-280.146 *pow(y,5.)+24.2578 *pow(y,6.)-0.876364 *pow(y,7.)+0.0108531 *pow(y,8.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(30.7098 -49.4572 *y+27.1065 *pow(y,2.)-5.81344 *pow(y,3.)+0.437775 *pow(y,4.)-0.00239131 *pow(y,5.)+3.39299e-6 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(37.4388 -61.4369 *y+35.8547 *pow(y,2.)-8.83604 *pow(y,3.)+0.837391 *pow(y,4.)-0.0138598 *pow(y,5.)+0.0000634027 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(34.5073 -56.5658 *y+32.6942 *pow(y,2.)-7.88295 *pow(y,3.)+0.726092 *pow(y,4.)-0.0118508 *pow(y,5.)+0.0000535522 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(39.2788 -77.5172 *y+54.4001 *pow(y,2.)-16.6524 *pow(y,3.)+2.25976 *pow(y,4.)-0.11185 *pow(y,5.)+0.000305286 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(55.7952 -109.485 *y+77.88 *pow(y,2.)-24.8458 *pow(y,3.)+3.59017 *pow(y,4.)-0.195215 *pow(y,5.)+0.00157563 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(59.3717 -115.437 *y+81.0835 *pow(y,2.)-25.3946 *pow(y,3.)+3.59089 *pow(y,4.)-0.191343 *pow(y,5.)+0.00152474 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(39.2788 -77.5172 *y+54.4001 *pow(y,2.)-16.6524 *pow(y,3.)+2.25976 *pow(y,4.)-0.11185 *pow(y,5.)+0.000305286 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(55.7952 -109.485 *y+77.88 *pow(y,2.)-24.8458 *pow(y,3.)+3.59017 *pow(y,4.)-0.195215 *pow(y,5.)+0.00157563 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(59.3717 -115.437 *y+81.0835 *pow(y,2.)-25.3946 *pow(y,3.)+3.59089 *pow(y,4.)-0.191343 *pow(y,5.)+0.00152474 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(50.2389 -117.386 *y+101.436 *pow(y,2.)-40.7468 *pow(y,3.)+8.12019 *pow(y,4.)-0.768511 *pow(y,5.)+0.0274682 *pow(y,6.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(83.1516 -189.881 *y+161.079 *pow(y,2.)-63.7554 *pow(y,3.)+12.4636 *pow(y,4.)-1.14312 *pow(y,5.)+0.039156 *pow(y,6.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(102.152 -229.779 *y+191.739 *pow(y,2.)-74.5327 *pow(y,3.)+14.3044 *pow(y,4.)-1.28895 *pow(y,5.)+0.0434125 *pow(y,6.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
            return np.exp(-2.*y)*(-177.322+386.798 *y-300.422 *pow(y,2.)+102.903 *pow(y,3.)-15.7025 *pow(y,4.)+0.88637 *pow(y,5.)-0.00302303 *pow(y,6.)+1.69649e-6 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-195.79+430.809 *y-341.448 *pow(y,2.)+122.197 *pow(y,3.)-20.0864 *pow(y,4.)+1.28975 *pow(y,5.)-0.0139876 *pow(y,6.)+0.0000317013 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-187.968+413.358 *y-326.309 *pow(y,2.)+115.653 *pow(y,3.)-18.7302 *pow(y,4.)+1.18206 *pow(y,5.)-0.0124973 *pow(y,6.)+0.0000267761 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0
    
    def FMPhi2pn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-226.801+577.066 *y-545.176 *pow(y,2.)+240.993 *pow(y,3.)-53.0006 *pow(y,4.)+5.54739 *pow(y,5.)-0.220344 *pow(y,6.)+0.000152643 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-291.786+735.781 *y-689.14 *pow(y,2.)+302.346 *pow(y,3.)-65.9297 *pow(y,4.)+6.83158 *pow(y,5.)-0.272999 *pow(y,6.)+0.000787813 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-323.409+809.866 *y-751.929 *pow(y,2.)+326.415 *pow(y,3.)-70.3157 *pow(y,4.)+7.18854 *pow(y,5.)-0.282787 *pow(y,6.)+0.000762371 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2np(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-243.719+590.938 *y-516.719 *pow(y,2.)+205.578 *pow(y,3.)-38.3324 *pow(y,4.)+3.01389 *pow(y,5.)-0.0611324 *pow(y,6.)+0.000152643 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-269.17+660.657 *y-589.938 *pow(y,2.)+244.605 *pow(y,3.)-48.9524 *pow(y,4.)+4.32244 *pow(y,5.)-0.118838 *pow(y,6.)+0.000787813 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-258.423+635.264 *y-567.834 *pow(y,2.)+235.143 *pow(y,3.)-46.91 *pow(y,4.)+4.13791 *pow(y,5.)-0.115134 *pow(y,6.)+0.000762371 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
            return 0

    def FMPhi2nn(self, ER):
        y = self.Y(ER)
        if self.shell_model=="Fitz":
            return np.exp(-2.*y)*(-311.725+868.999 *y-909.624 *pow(y,2.)+455.227 *pow(y,3.)-116.673 *pow(y,4.)+15.0231 *pow(y,5.)-0.852789 *pow(y,6.)+0.0137341 *pow(y,7.))
        elif self.shell_model=="JUN45":
            return np.exp(-2.*y)*(-401.145+1113.46 *y-1159.78 *pow(y,2.)+578.46 *pow(y,3.)-147.99 *pow(y,4.)+19.0734 *pow(y,5.)-1.09917 *pow(y,6.)+0.019578 *pow(y,7.))
        elif self.shell_model=="jj44b":
            return np.exp(-2.*y)*(-444.63+1228.65 *y-1274.04 *pow(y,2.)+632.912 *pow(y,3.)-161.371 *pow(y,4.)+20.761 *pow(y,5.)-1.1986 *pow(y,6.)+0.0217063 *pow(y,7.))
        else:
            print(self.shell_model+" doesn't exist for Ge.")
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
class Ge(Target):

    def __init__(self, shell_model="Fitz"):
        # Choose the shell model on initialisation to make it easier to customise on the detector side rather than when writing models.
        # Initialise all the isotopes: 70Ge (20%), 72Ge (27.4%), 73Ge (7.76%), 74Ge (36.5%), 76Ge (7.75%).
        # We'll assume we always want to use the same shell model, but somewhat trivial to tweak if we decide to change that
        self.ge70 = Ge70(shell_model)
        self.ge72 = Ge72(shell_model)
        self.ge73 = Ge73(shell_model)
        self.ge74 = Ge74(shell_model)
        self.ge76 = Ge76(shell_model)
        
    def A(self):
        return 0.205*self.ge70.A()+0.274*self.ge72.A()+0.0776*self.ge73.A()+0.365*self.ge74.A()+0.0775*self.ge76.A()

    def FMpp(self, ER):
        return 0.205*self.ge70.FMpp(ER)+0.274*self.ge72.FMpp(ER)+0.0776*self.ge73.FMpp(ER)+0.365*self.ge74.FMpp(ER)+0.0775*self.ge76.FMpp(ER)

    def FMpn(self, ER):
        return 0.205*self.ge70.FMpn(ER)+0.274*self.ge72.FMpn(ER)+0.0776*self.ge73.FMpn(ER)+0.365*self.ge74.FMpn(ER)+0.0775*self.ge76.FMpn(ER)
        
    def FMnp(self, ER):
        return 0.205*self.ge70.FMnp(ER)+0.274*self.ge72.FMnp(ER)+0.0776*self.ge73.FMnp(ER)+0.365*self.ge74.FMnp(ER)+0.0775*self.ge76.FMnp(ER)
        
    def FMnn(self, ER):
        return 0.205*self.ge70.FMnn(ER)+0.274*self.ge72.FMnn(ER)+0.0776*self.ge73.FMnn(ER)+0.365*self.ge74.FMnn(ER)+0.0775*self.ge76.FMnn(ER)

    def FS1pp(self, ER):
        return 0.205*self.ge70.FS1pp(ER)+0.274*self.ge72.FS1pp(ER)+0.0776*self.ge73.FS1pp(ER)+0.365*self.ge74.FS1pp(ER)+0.0775*self.ge76.FS1pp(ER)
    
    def FS1pn(self, ER):
        return 0.205*self.ge70.FS1pn(ER)+0.274*self.ge72.FS1pn(ER)+0.0776*self.ge73.FS1pn(ER)+0.365*self.ge74.FS1pn(ER)+0.0775*self.ge76.FS1pn(ER)

    def FS1np(self, ER):
        return 0.205*self.ge70.FS1np(ER)+0.274*self.ge72.FS1np(ER)+0.0776*self.ge73.FS1np(ER)+0.365*self.ge74.FS1np(ER)+0.0775*self.ge76.FS1np(ER)

    def FS1nn(self, ER):
        return 0.205*self.ge70.FS1nn(ER)+0.274*self.ge72.FS1nn(ER)+0.0776*self.ge73.FS1nn(ER)+0.365*self.ge74.FS1nn(ER)+0.0775*self.ge76.FS1nn(ER)

    def FS2pp(self, ER):
        return 0.205*self.ge70.FS2pp(ER)+0.274*self.ge72.FS2pp(ER)+0.0776*self.ge73.FS2pp(ER)+0.365*self.ge74.FS2pp(ER)+0.0775*self.ge76.FS2pp(ER)
    
    def FS2pn(self, ER):
        return 0.205*self.ge70.FS2pn(ER)+0.274*self.ge72.FS2pn(ER)+0.0776*self.ge73.FS2pn(ER)+0.365*self.ge74.FS2pn(ER)+0.0775*self.ge76.FS2pn(ER)

    def FS2np(self, ER):
        return 0.205*self.ge70.FS2np(ER)+0.274*self.ge72.FS2np(ER)+0.0776*self.ge73.FS2np(ER)+0.365*self.ge74.FS2np(ER)+0.0775*self.ge76.FS2np(ER)

    def FS2nn(self, ER):
        return 0.205*self.ge70.FS2nn(ER)+0.274*self.ge72.FS2nn(ER)+0.0776*self.ge73.FS2nn(ER)+0.365*self.ge74.FS2nn(ER)+0.0775*self.ge76.FS2nn(ER)

    def FPhi2pp(self, ER):
        return 0.205*self.ge70.FPhi2pp(ER)+0.274*self.ge72.FPhi2pp(ER)+0.0776*self.ge73.FPhi2pp(ER)+0.365*self.ge74.FPhi2pp(ER)+0.0775*self.ge76.FPhi2pp(ER)
    
    def FPhi2pn(self, ER):
        return 0.205*self.ge70.FPhi2pn(ER)+0.274*self.ge72.FPhi2pn(ER)+0.0776*self.ge73.FPhi2pn(ER)+0.365*self.ge74.FPhi2pn(ER)+0.0775*self.ge76.FPhi2pn(ER)

    def FPhi2np(self, ER):
        return 0.205*self.ge70.FPhi2np(ER)+0.274*self.ge72.FPhi2np(ER)+0.0776*self.ge73.FPhi2np(ER)+0.365*self.ge74.FPhi2np(ER)+0.0775*self.ge76.FPhi2np(ER)

    def FPhi2nn(self, ER):
        return 0.205*self.ge70.FPhi2nn(ER)+0.274*self.ge72.FPhi2nn(ER)+0.0776*self.ge73.FPhi2nn(ER)+0.365*self.ge74.FPhi2nn(ER)+0.0775*self.ge76.FPhi2nn(ER)
        
    def FPhipp(self, ER):
        return 0.205*self.ge70.FPhipp(ER)+0.274*self.ge72.FPhipp(ER)+0.0776*self.ge73.FPhipp(ER)+0.365*self.ge74.FPhipp(ER)+0.0775*self.ge76.FPhipp(ER)
    
    def FPhipn(self, ER):
        return 0.205*self.ge70.FPhipn(ER)+0.274*self.ge72.FPhipn(ER)+0.0776*self.ge73.FPhipn(ER)+0.365*self.ge74.FPhipn(ER)+0.0775*self.ge76.FPhipn(ER)

    def FPhinp(self, ER):
        return 0.205*self.ge70.FPhinp(ER)+0.274*self.ge72.FPhinp(ER)+0.0776*self.ge73.FPhinp(ER)+0.365*self.ge74.FPhinp(ER)+0.0775*self.ge76.FPhinp(ER)

    def FPhinn(self, ER):
        return 0.205*self.ge70.FPhinn(ER)+0.274*self.ge72.FPhinn(ER)+0.0776*self.ge73.FPhinn(ER)+0.365*self.ge74.FPhinn(ER)+0.0775*self.ge76.FPhinn(ER)

    def FDpp(self, ER):
        return 0.205*self.ge70.FDpp(ER)+0.274*self.ge72.FDpp(ER)+0.0776*self.ge73.FDpp(ER)+0.365*self.ge74.FDpp(ER)+0.0775*self.ge76.FDpp(ER)
    
    def FDpn(self, ER):
        return 0.205*self.ge70.FDpn(ER)+0.274*self.ge72.FDpn(ER)+0.0776*self.ge73.FDpn(ER)+0.365*self.ge74.FDpn(ER)+0.0775*self.ge76.FDpn(ER)

    def FDnp(self, ER):
        return 0.205*self.ge70.FDnp(ER)+0.274*self.ge72.FDnp(ER)+0.0776*self.ge73.FDnp(ER)+0.365*self.ge74.FDnp(ER)+0.0775*self.ge76.FDnp(ER)

    def FDnn(self, ER):
        return 0.205*self.ge70.FDnn(ER)+0.274*self.ge72.FDnn(ER)+0.0776*self.ge73.FDnn(ER)+0.365*self.ge74.FDnn(ER)+0.0775*self.ge76.FDnn(ER)

    def FMPhi2pp(self, ER):
        return 0.205*self.ge70.FMPhi2pp(ER)+0.274*self.ge72.FMPhi2pp(ER)+0.0776*self.ge73.FMPhi2pp(ER)+0.365*self.ge74.FMPhi2pp(ER)+0.0775*self.ge76.FMPhi2pp(ER)
    
    def FMPhi2pn(self, ER):
        return 0.205*self.ge70.FMPhi2pn(ER)+0.274*self.ge72.FMPhi2pn(ER)+0.0776*self.ge73.FMPhi2pn(ER)+0.365*self.ge74.FMPhi2pn(ER)+0.0775*self.ge76.FMPhi2pn(ER)

    def FMPhi2np(self, ER):
        return 0.205*self.ge70.FMPhi2np(ER)+0.274*self.ge72.FMPhi2np(ER)+0.0776*self.ge73.FMPhi2np(ER)+0.365*self.ge74.FMPhi2np(ER)+0.0775*self.ge76.FMPhi2np(ER)

    def FMPhi2nn(self, ER):
        return 0.205*self.ge70.FMPhi2nn(ER)+0.274*self.ge72.FMPhi2nn(ER)+0.0776*self.ge73.FMPhi2nn(ER)+0.365*self.ge74.FMPhi2nn(ER)+0.0775*self.ge76.FMPhi2nn(ER)

    def FS1Dpp(self, ER):
        return 0.205*self.ge70.FS1Dpp(ER)+0.274*self.ge72.FS1Dpp(ER)+0.0776*self.ge73.FS1Dpp(ER)+0.365*self.ge74.FS1Dpp(ER)+0.0775*self.ge76.FS1Dpp(ER)

    def FS1Dpn(self, ER):
        return 0.205*self.ge70.FS1Dpn(ER)+0.274*self.ge72.FS1Dpn(ER)+0.0776*self.ge73.FS1Dpn(ER)+0.365*self.ge74.FS1Dpn(ER)+0.0775*self.ge76.FS1Dpn(ER)

    def FS1Dnp(self, ER):
        return 0.205*self.ge70.FS1Dnp(ER)+0.274*self.ge72.FS1Dnp(ER)+0.0776*self.ge73.FS1Dnp(ER)+0.365*self.ge74.FS1Dnp(ER)+0.0775*self.ge76.FS1Dnp(ER)

    def FS1Dnn(self, ER):
        return 0.205*self.ge70.FS1Dnn(ER)+0.274*self.ge72.FS1Dnn(ER)+0.0776*self.ge73.FS1Dnn(ER)+0.365*self.ge74.FS1Dnn(ER)+0.0775*self.ge76.FS1Dnn(ER)  

## Will define these three later...
    def eNL(self):
        pass

    def eTransE_E(self):
        ## This is going to be potentially annoying to define with multi-isotopes and the way its being done in the back end... need to think about how to deal with ith
        pass
    
    def sigma_PE(self,E_gam):
        pass
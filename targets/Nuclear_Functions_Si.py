
## Nuclear Functions For Si

## Wrote the functions with the second argument denoting the shell model interaction used. F= Fitzpatrick ; the rest are based on the shell model names 



 class Si28(Target):
    def A(self):
    return 28

## Fitzpatrick (USD)

# M channel
        def FMpp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))
        
        def FMpn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))

        def FMnp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))

        def FMnn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(195.992 -335.983*y+196.671 *pow(y,2.)-45.1536 *pow(y,3.)+3.53988 *pow(y,4.))

# Sigma'' (Sigma2)

        def FS2pp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FS2pn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS2np(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS2nn(self, F, ER):
           y = self.Y(ER)
            return 0

# Sigma' (Sigma1)

        def FS1pp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FS1pn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1np(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1nn(self, F, ER):
           y = self.Y(ER)
            return 0


#Phi'' (Phi2)
        def FPhi2pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))
        
        def FPhi2pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))

        def FPhi2np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))

        def FPhi2nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(5.8049 -4.64392*y+0.928784 *pow(y,2.))
        
#Phi' (Phi)
        def FPhipp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FPhipn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FPhinp(self, F, ER):
           y = self.Y(ER)
            return 0

        def FPhinn(self, F, ER):
           y = self.Y(ER)
            return 0

# Delta (D)
        def FDpp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FDpn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FDnp(self, F, ER):
           y = self.Y(ER)
            return 0

        def FDnn(self, F, ER):
           y = self.Y(ER)
            return 0


# M Phi2
        def FMPhi2pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))
        
        def FMPhi2pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))

        def FMPhi2np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))

        def FMPhi2nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-33.7301+42.4032*y-16.0975 *pow(y,2.)+1.81323 *pow(y,3.))

# Sigma1 Delta
        def FS1Dpp(self, F, ER):
           y = self.Y(ER)
            return 0
    
        def FS1Dpn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnp(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnn(self, F, ER):
           y = self.Y(ER)
            return 0   

        
##USDB

# M channel

    def FMpp(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))

    def FMpn(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))

    def FMnp(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))

    def FMnn(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(195.998 -335.997*y+196.8 *pow(y,2.)-45.2579 *pow(y,3.)+3.55606 *pow(y,4.))
    

# Sigma'' (Sigma2)

        def FS2pp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FS2pn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS2np(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS2nn(self, USDB, ER):
           y = self.Y(ER)
            return 0

# Sigma' (Sigma1)

        def FS1pp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FS1pn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1np(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1nn(self, USDB, ER):
           y = self.Y(ER)
            return 0

#Phi'' (Phi2)
        def FPhi2pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))
        
        def FPhi2pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))

        def FPhi2np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))

        def FPhi2nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(6.14601 -4.91684*y+0.983372 *pow(y,2.))

#Phi' (Phi)
        def FPhipp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FPhipn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FPhinp(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FPhinn(self, USDB, ER):
           y = self.Y(ER)
            return 0

# Delta (D)
        def FDpp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FDpn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FDnp(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FDnn(self, USDB, ER):
           y = self.Y(ER)
            return 0


# M Phi2
        def FMPhi2pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))
        
        def FMPhi2pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))

        def FMPhi2np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))

        def FMPhi2nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-34.7078+43.6326*y-16.5748 *pow(y,2.)+1.87001 *pow(y,3.))

# Sigma1 Delta
        def FS1Dpp(self, USDB, ER):
           y = self.Y(ER)
            return 0
    
        def FS1Dpn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnp(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnn(self, USDB, ER):
           y = self.Y(ER)
            return 0   



class Si29(Target):
    def A(self):
        return 29

## Fitzpatrick (USD)
    def FMpp(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(195.994 -335.987*y+195.357*pow(y,2.)-44.0261 *pow(y,3.)+3.36526 *pow(y,4.))

    def FMpn(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(209.992 -366.649*y+219.759 *pow(y,2.)-52.1025 *pow(y,3.)+4.22607 *pow(y,4.))

    def FMnp(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(209.992 -366.649*y+219.759 *pow(y,2.)-52.1025 *pow(y,3.)+4.22607 *pow(y,4.))

    def FMnn(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(224.99 -399.978*y+246.876 *pow(y,2.)-61.4301 *pow(y,3.)+5.30706 *pow(y,4.))


# Sigma'' (Sigma2)

        def FS2pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.0000165528 +0.0000306321*y-0.0000519747 *pow(y,2.)-0.0000612043 *pow(y,3.)+0.0000660817 *pow(y,4.))
        
        def FS2pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.00108185-0.000980401*y+0.00134241 *pow(y,2.)-0.000816813 *pow(y,3.)+0.00167485 *pow(y,4.))

        def FS2np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.00108185-0.000980401*y+0.00134241 *pow(y,2.)-0.000816813 *pow(y,3.)+0.00167485 *pow(y,4.))

        def FS2nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.0707069 -0.0026951*y+0.109597 *pow(y,2.)-0.00208824 *pow(y,3.)+0.0424495 *pow(y,4.))


# Sigma' (Sigma1)

        def FS1pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)* (0.0000331056 -0.000163055*y+0.000199284 *pow(y,2.)+3.66554e-6 *pow(y,3.)+1.67306e-8 *pow(y,4.))
        
        def FS1pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.0021637+0.0096352*y-0.0133306 *pow(y,2.)+0.00673258 *pow(y,3.)+0.0000623427 *pow(y,4.))

        def FS1np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.0021637+0.0096352*y-0.0133306 *pow(y,2.)+0.00673258 *pow(y,3.)+0.0000623427 *pow(y,4.))

        def FS1nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.141414 -0.56296*y+0.922776 *pow(y,2.)-0.721542 *pow(y,3.)+0.232306 *pow(y,4.))


#Phi'' (Phi2)
        def FPhi2pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(7.52516 -6.02013*y+1.20403 *pow(y,2.))
        
        def FPhi2pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(7.62376 -6.09901*y+1.2198 *pow(y,2.))

        def FPhi2np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(7.62376 -6.09901*y+1.2198 *pow(y,2.))

        def FPhi2nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(7.72367 -6.17894*y+1.23579 *pow(y,2.))


#Phi' (Phi)
        def FPhipp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FPhipn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FPhinp(self, F, ER):
           y = self.Y(ER)
            return 0

        def FPhinn(self, F, ER):
           y = self.Y(ER)
            return 0

# Delta (D)
        def FDpp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.000176249 -0.000140999*y+0.0000281999 *pow(y,2.))
        
        def FDpn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.00328811 -0.00263049*y+0.000526098 *pow(y,2.))

        def FDnp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.00328811 -0.00263049*y+0.000526098 *pow(y,2.))

        def FDnn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.0613433 -0.0490746*y+0.00981492 *pow(y,2.))

# M Phi2
        def FMPhi2pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-38.4042+48.2793*y-18.1994 *pow(y,2.)+2.01292 *pow(y,3.))
        
        def FMPhi2pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-38.9075+48.912*y-18.4378 *pow(y,2.)+2.0393 *pow(y,3.))

        def FMPhi2np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-41.1471+53.0337*y-20.9495 *pow(y,2.)+2.52781 *pow(y,3.))

        def FMPhi2nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-41.6863+53.7286*y-21.224 *pow(y,2.)+2.56094 *pow(y,3.))

# Sigma1 Delta
        def FS1Dpp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.0000763861 -0.000218666*y+0.0000735275 *pow(y,2.)+6.86877e-7 *pow(y,3.))
    
        def FS1Dpn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.00142506 -0.00407945*y+0.00137173 *pow(y,2.)+0.0000128144 *pow(y,3.))

        def FS1Dnp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.0049924+0.0119342*y-0.0103736 *pow(y,2.)+0.00255949 *pow(y,3.))

        def FS1Dnn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.0931385+0.222645*y-0.193531 *pow(y,2.)+0.04775 *pow(y,3.))



## USDB

    def FMpp(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(195.998 -335.997*y+195.406 *pow(y,2.)-44.0637 *pow(y,3.)+3.37088 *pow(y,4.))

    def FMpn(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(209.996 -366.658*y+219.796 *pow(y,2.)-52.1321 *pow(y,3.)+4.23088 *pow(y,4.))

    def FMnp(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(209.996 -366.658*y+219.796 *pow(y,2.)-52.1321 *pow(y,3.)+4.23088 *pow(y,4.))

    def FMnn(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(224.992 -399.984*y+246.9 *pow(y,2.)-61.4494 *pow(y,3.)+5.3103 *pow(y,4.))


# Sigma'' (Sigma2)

        def FS2pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.00102417 -0.00375031*y+0.00485202 *pow(y,2.)-0.00259768 *pow(y,3.)+0.000491371 *pow(y,4.))
        
        def FS2pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.00992846 -0.0206121*y+0.0187405 *pow(y,2.)-0.0152473 *pow(y,3.)+0.00513046 *pow(y,4.))

        def FS2np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.00992846 -0.0206121*y+0.0187405 *pow(y,2.)-0.0152473 *pow(y,3.)+0.00513046 *pow(y,4.))

        def FS2nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.0962483 -0.0471923*y+0.149393 *pow(y,2.)-0.0352067 *pow(y,3.)+0.0535677 *pow(y,4.))


# Sigma' (Sigma1)

        def FS1pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.00204833 -0.00444303*y+0.00359312 *pow(y,2.)-0.00128387 *pow(y,3.)+0.000171036 *pow(y,4.))
        
        def FS1pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.0198569 -0.0588155*y+0.0694612 *pow(y,2.)-0.0360334 *pow(y,3.)+0.00673044 *pow(y,4.))

        def FS1np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.0198569 -0.0588155*y+0.0694612 *pow(y,2.)-0.0360334 *pow(y,3.)+0.00673044 *pow(y,4.))

        def FS1nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(0.192497 -0.722794*y+1.13008 *pow(y,2.)-0.847818 *pow(y,3.)+0.26485 *pow(y,4.))

#Phi'' (Phi2)
        def FPhi2pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(8.0029 -6.40235*y+1.28048 *pow(y,2.))
        
        def FPhi2pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(8.17622 -6.541*y+1.30821 *pow(y,2.))

        def FPhi2np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(8.17622 -6.541*y+1.30821 *pow(y,2.))

        def FPhi2nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(8.35329 -6.68266*y+1.33654 *pow(y,2.))

#Phi' (Phi)
        def FPhipp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FPhipn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FPhinp(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FPhinn(self, USDB, ER):
           y = self.Y(ER)
            return 0

# Delta (D)
        def FDpp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*
        
        def FDpn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

        def FDnp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

        def FDnn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

# M Phi2
        def FMPhi2pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-39.6061+49.7905*y-18.7733 *pow(y,2.)+2.07764 *pow(y,3.))
        
        def FMPhi2pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-40.4632+50.868*y-19.1796 *pow(y,2.)+2.1226 *pow(y,3.))

        def FMPhi2np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-42.4351+54.6941*y-21.6074 *pow(y,2.)+2.60776 *pow(y,3.))

        def FMPhi2nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-43.3534+55.8777*y-22.075 *pow(y,2.)+2.6642 *pow(y,3.))

# Sigma1 Delta
        def FS1Dpp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.0000504142+0.000074823*y-0.000036418 *pow(y,2.)+5.82204e-6 *pow(y,3.))
    
        def FS1Dpn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.0104582+0.0155218*y-0.00755479 *pow(y,2.)+0.00120776 *pow(y,3.))

        def FS1Dnp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.000489559+0.00111492*y-0.0009418 *pow(y,2.)+0.000229665 *pow(y,3.))

        def FS1Dnn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-0.101557+0.231286*y-0.195373 *pow(y,2.)+0.0476431 *pow(y,3.))





class Si30(Target):
    def A(self):
        return 30

## Fitzpatrick (USD)
    
    def FMpp(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(195.996 -335.991*y+194.701 *pow(y,2.)-43.4623 *pow(y,3.)+3.27957 *pow(y,4.))

    def FMpn(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(223.993 -397.318*y+242.474 *pow(y,2.)-58.709 *pow(y,3.)+4.8518 *pow(y,4.))

    def FMnp(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(223.993 -397.318*y+242.474 *pow(y,2.)-58.709 *pow(y,3.)+4.8518 *pow(y,4.))

    def FMnn(self, F, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(255.99 -469.311*y+300.829 *pow(y,2.)-78.5857 *pow(y,3.)+7.17776 *pow(y,4.))


# Sigma'' (Sigma2)

        def FS2pp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FS2pn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS2np(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS2nn(self, F, ER):
           y = self.Y(ER)
            return 0

# Sigma' (Sigma1)

        def FS1pp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FS1pn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1np(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1nn(self, F, ER):
           y = self.Y(ER)
            return 0

#Phi'' (Phi2)
        def FPhi2pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*
        
        def FPhi2pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

        def FPhi2np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

        def FPhi2nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

#Phi' (Phi)
        def FPhipp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FPhipn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FPhinp(self, F, ER):
           y = self.Y(ER)
            return 0

        def FPhinn(self, F, ER):
           y = self.Y(ER)
            return 0

# Delta (D)
        def FDpp(self, F, ER):
           y = self.Y(ER)
            return 0
        
        def FDpn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FDnp(self, F, ER):
           y = self.Y(ER)
            return 0

        def FDnn(self, F, ER):
           y = self.Y(ER)
            return 0


# M Phi2
        def FMPhi2pp(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-41.0858+51.6505*y-19.4011 *pow(y,2.)+2.12587 *pow(y,3.))
        
        def FMPhi2pn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-33.7662+42.4487*y-15.9447 *pow(y,2.)+1.74713 *pow(y,3.))

        def FMPhi2np(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-46.9548+61.8234*y-25.0791 *pow(y,2.)+3.14501 *pow(y,3.))

        def FMPhi2nn(self, F, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-38.5896+50.8093*y-20.6112 *pow(y,2.)+2.58471 *pow(y,3.))

# Sigma1 Delta
        def FS1Dpp(self, F, ER):
           y = self.Y(ER)
            return 0
    
        def FS1Dpn(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnp(self, F, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnn(self, F, ER):
           y = self.Y(ER)
            return 0   





## USDB

#  M 

    def FMpp(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(196.001 -336.007*y+194.784 *pow(y,2.)-43.5249 *pow(y,3.)+3.2888 *pow(y,4.))

    def FMpn(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(223.995 -397.326*y+242.532 *pow(y,2.)-58.7566 *pow(y,3.)+4.86 *pow(y,4.))

    def FMnp(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(223.995 -397.326*y+242.532 *pow(y,2.)-58.7566 *pow(y,3.)+4.86 *pow(y,4.))

    def FMnn(self, USDB, ER):
       y = self.Y(ER)
        return np.exp(-2.*y)*(255.986 -469.307*y+300.853 *pow(y,2.)-78.608 *pow(y,3.)+7.18185 *pow(y,4.))


# Sigma'' (Sigma2)

        def FS2pp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FS2pn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS2np(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS2nn(self, USDB, ER):
           y = self.Y(ER)
            return 0

# Sigma' (Sigma1)

        def FS1pp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FS1pn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1np(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1nn(self, USDB, ER):
           y = self.Y(ER)
            return 0


#Phi'' (Phi2)
        def FPhi2pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*
        
        def FPhi2pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

        def FPhi2np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

        def FPhi2nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*

#Phi' (Phi)
        def FPhipp(self, USDB, ER):
           y = self.Y(ER)
            return 0
        
        def FPhipn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FPhinp(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FPhinn(self, USDB, ER):
           y = self.Y(ER)
            return 0

# Delta (D)
        def FDpp(self, USDB, ER):
           y = self.Y(ER)
            return 0
    
        def FDpn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FDnp(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FDnn(self, USDB, ER):
           y = self.Y(ER)
            return 0



# M Phi2
        def FMPhi2pp(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-42.3041+53.1822*y-19.9837 *pow(y,2.)+2.19178 *pow(y,3.))
        
        def FMPhi2pn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-35.1452+44.1825*y-16.602 *pow(y,2.)+1.82088 *pow(y,3.))

        def FMPhi2np(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-48.3476+63.6577*y-25.8258 *pow(y,2.)+3.23935 *pow(y,3.))

        def FMPhi2nn(self, USDB, ER):
           y = self.Y(ER)
            return np.exp(-2.*y)*(-40.166+52.8853*y-21.4555 *pow(y,2.)+2.69117 *pow(y,3.))


# Sigma1 Delta
        def FS1Dpp(self, USDB, ER):
           y = self.Y(ER)
            return 0
    
        def FS1Dpn(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnp(self, USDB, ER):
           y = self.Y(ER)
            return 0

        def FS1Dnn(self, USDB, ER):
           y = self.Y(ER)
            return 0   


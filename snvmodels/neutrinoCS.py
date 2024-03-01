import numpy as np
from constants import *

def Sig(self,Target,Ev,ER):
    """
    Differential cross section as a function of recoil energy ER in [eV]
    Target: target nucleus
    mX: DM mass [eV]
    ER: recoil energy [eV]
    Sig: cross section [cm2/eV]
    """
    ##inverse beta decay cross section (IBD) https://arxiv.org/pdf/astro-ph/0302055.pdf
    costhetaC = 0.9746
    xi=3.706
    MA2=1.014*GeV*GeV
    #epsilon = Ev/mp
    t= (mn**2-mproton**2-2*mproton*(Ev-ER))
    f1=1+2.41*t/GeV**2
    f2=xi*(1+3.21*t/GeV**2)
    g0=-1.270
    g1=g0/(1-t/MA2)**2
    A = M_nucleus**2*(f1**2-g1**2)*(t-me**2)-M_nucleus**2*Delta_np**2*(f1**2+g1**2)-2*me**2*M_nucleus*Delta_np*g1*(f1+f2)
    B=t*g1*(f1+f2)
    C=(f1**2+g1**2)/4
    M2 = A-(2*mproton*(Ev+ER)-me**2)*B+(2*mproton*(Ev+ER)-me**2)**2*C
    #print("M2",M2)
    cross_sec = Gf**2*costhetaC**2/(2*np.pi*(2*mproton*Ev)**2)*M2*0.389379*1e-27/1e9
    return cross_sec
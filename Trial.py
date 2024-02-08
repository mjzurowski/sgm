import numpy as np
import matplotlib.pyplot as plt
import BSWDet
#import BSW
from constants import *

Det = BSWDet.DAMA_BSW
mX = 690*keV
sig = 1
#Model = BSW.SIWIMPBSW()

E = np.arange(10,100,0.1) # observed energy, units of keV_ee
plt.plot(E,[Det.dRdE(e,mX=mX,sig=sig)for e in E])
plt.ylabel(r"Observation rate [cpd/kg/keV$_{\rm ee}$]")
plt.xlabel(r"Energy [keV$_{\rm ee}$]")
#plt.xlim(0,16)
#plt.ylim(0,0.026)
plt.grid()
plt.show()
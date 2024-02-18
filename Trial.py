import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from models.wimp import SIWIMP
from models.wimp import SIWIMP_Helm
from models.inelastic import SIInel
from detectors.dama import DAMA
from detectors.cosine import COSINE
from veldists import VelDist
from constants import *
import BSW2
from BSWDet import  DAMA_BSW
## Step 1
Det = DAMA_BSW()
Model1 = SIWIMP()
Model2 = SIInel()
Dist = VelDist("modSHM",0.3) # standard halo model
Dist2 = VelDist("modShards",0.3) # standard halo model + shards
Model = BSW2.SIWIMPBSW2()

E = np.arange(10,1000,0.1) # observed energy, units of keV_ee
plt.plot(E,[Det.dRdE(e, Model.dRdER,mX=690*keV,sig=1,VelDist=Dist) for e in E])
plt.ylabel(r"Observation rate [cpd/kg/keV$_{\rm ee}$]")
plt.xlabel(r"Energy [keV$_{\rm ee}$]")
#plt.xlim(0,16)
#plt.ylim(0,0.026)
plt.grid()
plt.show()
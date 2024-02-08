import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import BSW 
from veldists import VelDist
from detectors.dama import DAMA
from constants import *




ER = np.arange(10,100,0.1) # observed energy, units of keV_ee
plt.plot(ER,[Det.dRdE(DMRate,e)for e in ER])
plt.ylabel(r"Observation rate [cpd/kg/keV$_{\rm ee}$]")
plt.xlabel(r"Energy [keV$_{\rm ee}$]")
#plt.xlim(0,16)
#plt.ylim(0,0.026)
plt.grid()
plt.show()
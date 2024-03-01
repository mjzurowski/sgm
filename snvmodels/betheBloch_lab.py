import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.constants import c
import scienceplots

#Constantes
pi = 3.1415926535
Na = 6.022e23 	# mol^-1
re = 2.8179e-15 # m
me = 0.511 		#MeV
C0 = 0.1535 	#MeV g^-1 cm^2
density_lab = 0.863 #g/cm3
mu =  931.4941 # unified atomic mass

beta2 = lambda bg: bg**2/(bg**2 + 1) 
gamma = lambda bg: np.sqrt(1 + bg**2)

def C(Z,hvp):
	return -2*np.log(I(Z)/hvp) - 1

def C_shell(bg,Z):
	aux1 = (0.422377*bg**(-2) + 0.0304043*bg**(-4) - 0.00038106*bg**(-6))*1e-6*I(Z)**2
	aux2 = (3.850190*bg**(-2) - 0.1667989*bg**(-4) + 0.00157955*bg**(-6))*1e-9*I(Z)**3

	return aux1 + aux2

def I(Z):
	if Z<13:
		aux =  12*Z + 7
	else:
		aux = 9.76*Z + 58.8*Z**(-0.19)

	return aux*1e-6 
		
def Wmax(bg,M):
	return (2*me*bg**2)/(1 + 2*gamma(bg)*me/M + (me/M)**2)

def delta(bg,Z,M,delta0,x0,x1,a,m,hvp): #delta0,x0,x1,a,k 
	x = np.log10(bg)

	if x >= x1:
		aux = 4.6052*x + C(Z,hvp)
	elif x >= x0:
		aux = 4.6052*x + C(Z,hvp) + a*(x1 - x)**m
	else:
		aux = delta0*10**(2*(x-x0))

	return aux
	
	

def dEdx(bg,z,M,a,m,x0,x1,delta0,Z,A,hvp):
	aux = C0*(Z/A)*(z**2/beta2(bg))*(np.log(2*me*bg*bg*Wmax(bg,M)/I(Z)**2) - 2*beta2(bg) - delta(bg,Z,M,delta0,x0,x1,a,m,hvp)
  - C_shell(bg,Z)/Z)
	
	return aux

 
def plotBethe(ee,z,M):
	yy_mat=0
	#print("z=",z)
	data = np.genfromtxt('snvmodels/data_BetheBloch.txt',dtype=None,encoding=None)
	for ind in range(0,2):
		a = float(data[ind][1])
		m = float(data[ind][2])
		x0 = float(data[ind][3])
		x1 = float(data[ind][4])
		delta0 = float(data[ind][5])
		Z = float(data[ind][6])
		A = float(data[ind][7])
		hvp = float(data[ind][8])*1e-6
		rho = float(data[ind][9])
		weight = float(data[ind][10])
		xx = np.sqrt((ee/M)**2+2*ee/M)
		yy=weight/rho*dEdx(xx,z,M,a,m,x0,x1,delta0,Z,A,hvp)*rho
		yy_mat+=yy
	return yy_mat*density_lab
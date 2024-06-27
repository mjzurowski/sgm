#Coefficent mappings between high energy and non-relativistic
#everything written in eV

import math
import sympy as sp
import numpy as np
from sympy import symbols

#ci(N, Coeff, mchi, Lambda) refer to the NR coefficients

#N=p,n mchi is DM mass Lambda is high energy physics coefficient Coeff is the name of the list which turns on or off certain cq and cg coefficients - to be set by user

##########################################
#Define vector/list which turns on or turns off different cq and cg coefficient
#bases on what we want. "1" turn on, "0" turns off

#Example, this turns on all of them
#CoeffList=[cq1,cq2,cq3,cq4,cq5,cq6,cq7,cq8,cq9,cq10,cg1,cg2,cg3,cg4]
#order of cq and cg terms needs to be like the above
#AllExample=[1,1,1,1,1,1,1,1,1,1,1,1,1,1]

#the name of the vector itself, i.e. AllExample, is what needs to be inputted in place of Coeff

#To turn on only one of them, e.g. cq1, do:
#Example=[1,0,0,0,0,0,0,0,0,0,0,0,0,0]
##########################################

#Define the cq and cg coefficients

Lambda =symbols('Lambda')

def cq1(q,Lambda):
    return m(q)/(Lambda**3)

def cq2(q,Lambda):
    return m(q)/(Lambda**3)

def cq3(q,Lambda):
    return m(q)/(Lambda**3)

def cq4(q,Lambda):
    return m(q)/(Lambda**3)

def cq5(Lambda):
    return 1/(Lambda**2)

def cq6(Lambda):
    return 1/(Lambda**2)

def cq7(Lambda):
    return 1/(Lambda**2)

def cq8(Lambda):
    return 1/(Lambda**2)

def cq9(Lambda):
    return 1/(Lambda**2)

def cq10(Lambda):
    return 1/(Lambda**2)

def cg1(Lambda):
    return None

def cg2(Lambda):
    return None

def cg3(Lambda):
    return None

def cg4(Lambda):
    return None


#Convert relativistic to quark-gloun operators
#quark and gluon operators are given by cqi and cgj
u,d,s,q,N,p,n,b,t,c =symbols('u d s q N p n b t c')

def fT(N,q):
    if N == p:
        if q==u:
            return 0.023
        if q==d:
            return 0.034
        if q==s:
            return 0.14
    if N==n:
        if q==u:
            return 0.019
        if q==d:
            return 0.041
        if q==s:
            return 0.14

def fTG(N):
    if N == p:
        return 1 - (fT(p,d) + fT(p,s) + fT(p,u))
    if N==n:
        return 1 - (fT(n,d) + fT(n,s) + fT(n,u))


def CDelta(N,q):
    if N == p:
        if q==u:
            return 0.77
        if q==d:
            return - 0.40
        if q==s:
            return - 0.12
        if q in [b,t,c]:
            return 0
    if N==n:
        if q==u:
            return - 0.40
        if q==d:
            return 0.77
        if q==s:
            return - 0.12
        if q in [b,t,c]:
            return 0

def Delta(N,q):
    if N == p:
        if q==u:
            return 0.84
        if q==d:
            return -0.23
        if q==s:
            return -0.05
        if q in [b,t,c]:
            return 0
    if N==n:
        if q==u:
            return -0.23
        if q==d:
            return 0.84
        if q==s:
            return -0.05
        if q in [b,t,c]:
            return 0

def m(q):
    if q==u:
        return 0.0023* 10**9  #eV
    if q==d:
        return 0.0048* 10**9  #eV
    if q==s:
        return 0.095* 10**9 #eV
    if q==b:
        return 4.18* 10**9  #eV
    if q==c:
        return 1.275* 10**9 #eV
    if q==t:
        return 173.210* 10**9  #eV

mBar = (1/m(u) + 1/m(d) + 1/m(s))**(-1) 




def C3(Coeff):
    return sum( Coeff[2]*cq3(q,Lambda)*mBar/m(q) for q in [u,d,s,t,b,c] )
def C4(Coeff):
    return sum( Coeff[3]*cq4(q,Lambda)*mBar/m(q) for q in [u,d,s,t,b,c] )


def Rc1(N, Lambda,Coeff):
    return sum( Coeff[0]*cq1(q,Lambda)*(mp/m(q)) * fT(N,q) for q in [u,d,s] ) + (2/27)*fTG(N) * ( sum( (mp/m(q))*Coeff[0]*cq1(q,Lambda) for q in [c,b,t] )  - Coeff[10]*cg1(Lambda) *mp)

def Rc2(N,Lambda,Coeff):
    return sum( Coeff[1]*cq2(q,Lambda)*(mp/m(q)) * fT(N,q) for q in [u,d,s] ) + (2/27)*fTG(N) * ( sum( (mp/m(q))*Coeff[1]*cq2(q,Lambda) for q in [c,b,t] )  - Coeff[11]*cg2(Lambda) *mp)

def Rc3(N,Lambda,Coeff):
    return sum( (mp/m(q))*CDelta(N,q)*( (Coeff[2]*cq3(q,Lambda)-C3(Coeff)) + Coeff[12]*cg3(Lambda)*mBar )  for q in [u,d,s] )

def Rc4(N,Lambda,Coeff):
    return sum( (mp/m(q))*CDelta(N,q)*( (Coeff[3]*cq4(q,Lambda)-C4(Coeff) ) + Coeff[13]*cg4(Lambda)*mBar )  for q in [u,d,s] )

def Rc5(N,Lambda,Coeff):
    if N==p:
        #supposed to have 2*cu5(Lambda) + cd5(Lambda), but cu5=cd5=cq5=1/Lambda^2
        return 2*Coeff[4]*cq5(Lambda) + Coeff[4]*cq5(Lambda)
    if N==n:
        #supposed to have cu5(Lambda) + 2*cd5(Lambda), but cu5=cd5=cq5=1/Lambda^2
        return Coeff[4]*cq5(Lambda) + 2*Coeff[4]*cq5(Lambda)

def Rc6(N,Lambda,Coeff):
    if N==p:
        return 2*Coeff[5]*cq6(Lambda) + Coeff[5]*cq6(Lambda)
    if N==n:
        return Coeff[5]*cq6(Lambda) + 2*Coeff[5]*cq6(Lambda)

def Rc7(N,Lambda,Coeff):
    return sum( Coeff[6]*cq7(Lambda)*CDelta(N,q) for q in [u,d,s,t,b,c] )

def Rc8(N,Lambda,Coeff):
    return sum( Coeff[7]*cq8(Lambda)*CDelta(N,q) for q in [u,d,s,t,b,c] )

def Rc9(N,Lambda,Coeff):
    return sum( Coeff[8]*cq9(Lambda)*Delta(N,q) for q in [u,d,s,t,b,c] )

def Rc10(N,Lambda,Coeff):
    return sum( Coeff[9]*cq10(Lambda)*Delta(N,q) for q in [u,d,s,t,b,c] )    



#Relativistic to NR conversion
#Relativistic coefficients are denoted Rcpi, Rcpj

#check these - do they go up to O15?

def c1(N,Coeff,mchi,Lambda):
    return 4*mchi*mp*Rc1(N, Lambda,Coeff) + 4*mchi*mp*Rc5(N, Lambda,Coeff)

def c3(N,Coeff,mchi,Lambda):
    return None #not defined

def c4(N,Coeff,mchi,Lambda):
    return -16*mchi*mp*Rc8(N, Lambda,Coeff) + 32* mchi*mp*Rc9(N, Lambda,Coeff)

def c5(N,Coeff,mchi,Lambda):
    return None #not defined

def c6(N,Coeff,mchi,Lambda):
    return 4*Rc4(N, Lambda,Coeff)

def c7(N,Coeff,mchi,Lambda):
    return -8*mchi*mp*Rc7(N, Lambda,Coeff)

def c8(N,Coeff,mchi,Lambda):
    return 8*mchi*mp*Rc6(N, Lambda,Coeff)

def c9(N,Coeff,mchi,Lambda):
    return 8*mchi*Rc6(N, Lambda,Coeff) + 8*mp*Rc7(N, Lambda,Coeff)

def c10(N,Coeff,mchi,Lambda):
    return 4*mchi*Rc3(N, Lambda,Coeff)-8*mp*Rc10(N, Lambda,Coeff)

def c11(N,Coeff,mchi,Lambda):
    return -4*mp*Rc2(N, Lambda,Coeff) + 8*mchi*Rc10(N, Lambda,Coeff)

def c12(N,Coeff,mchi,Lambda):
    return -32*mchi*mp*Rc10(N, Lambda,Coeff)





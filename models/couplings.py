import numpy as np
from constants import *

###################################################################################
### Defining quark level couplings where desired. See arxiv 1307.5955 for details
###################################################################################

## Quark masses
m_q = {
    "u":0.0023*GeV,
    "d":0.0048*GeV,
    "s":0.095*GeV,
    "b":4.18*GeV,
    "c":1.275*GeV,
    "t":173.210*GeV
    }

## Quark-Nucleon Delta param
## See Table 4 in 1307.5955
Delta_p = {
    "u":0.77,
    "d":-0.40,
    "s":-0.12,
    "b":0,
    "c":0,
    "t":0
    }

Delta_n = {
    "u":-0.40,
    "d":0.77,
    "s":-0.12,
    "b":0,
    "c":0,
    "t":0
    }

Delta = {
    "n":Delta_n,
    "p":Delta_p
    }

## Quark-Nucleon fT param
## See Table 4 in 1307.5955
fT_p = {
    "u":0.023,
    "d":0.034,
    "s":0.14,
    "b":0,
    "c":0,
    "t":0
    }

fT_n = {
    "u":0.019,
    "d":0.041,
    "s":0.14,
    "b":0,
    "c":0,
    "t":0
    }

fT = {
    "n":fT_n,
    "p":fT_p
    }

fTG = {
    "n":1 - (fT_n["u"] + fT_n["s"] + fT_n["d"]),
    "p":1 - (fT_p["u"] + fT_p["s"] + fT_p["d"])
}

## Quark-Nucleon delta (tensor coupling) params
## See 1307.5955 for proton defs, with n defined assuming isospin symmetry
delta_p = {
    "u":0.84,
    "d":-0.23,
    "s":-0.05,
    "b":0,
    "c":0,
    "t":0
    }

delta_n = {
    "u":-0.23,
    "d":0.84,
    "s":-0.05,
    "b":0,
    "c":0,
    "t":0
    }

delta = {
    "n":delta_n,
    "p":delta_p
    }

## Define a default dictionary for the quark couplings
## Assume we take the value in c_q and multiply by mq/Lambda^3 for 1-4, and 1/Lambda^2 for all others

def c_q():
    c_q = {
        "c1":0,
        "c2":0,
        "c3":0,
        "c4":0,
        "c5":0,
        "c6":0,
        "c7":0,
        "c8":0,
        "c9":0,
        "c10":0,
    }
    return c_q

mBar = 1/(1/m_q["u"]+1/m_q["d"]+1/m_q["s"])

###################################################
### Relativistic couplings
###################################################

def c1_N(cq):
    ## Define coupling, with the 1/Lambda^3 pulled out
    ## Input: dictionary of couplings
    ## Output units: [eV]
    if cq["c1"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s"]:
            cp += cq["c1"]*mp*fT["p"][q]
            cn += cq["c1"]*mp*fT["n"][q]
        for q in ["b","c","t"]:
            cp += (2/27)*cq["c1"]*mp*fTG["p"]
            cn += (2/27)*cq["c1"]*mp*fTG["n"]
        return [cp,cn]
        
def c2_N(cq):
    ## Define coupling, with the 1/Lambda^3 pulled out
    ## Input: dictionary of couplings
    ## Output units: [eV]
    if cq["c2"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s"]:
            cp += cq["c2"]*mp*fT["p"][q]
            cn += cq["c2"]*mp*fT["n"][q]
        for q in ["b","c","t"]:
            cp += (2/27)*cq["c2"]*mp*fTG["p"]
            cn += (2/27)*cq["c2"]*mp*fTG["n"]
        return [cp,cn]
    
def c3_N(cq):
    ## Define coupling, with the 1/Lambda^3 pulled out
    ## Input: dictionary of couplings
    ## Output units: [eV]
    if cq["c3"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s"]:
            cp += (cq["c3"]*mp/m_q[q])*(m_q[q]-6*mBar)*Delta["p"][q] # factor of six is because I'm assuming we sum over all quarks
            cn += (cq["c3"]*mp/m_q[q])*(m_q[q]-6*mBar)*Delta["n"][q] # factor of six is because I'm assuming we sum over all quarks
        return [cp,cn]
    
def c4_N(cq):
    ## Define coupling, with the 1/Lambda^3 pulled out
    ## Input: dictionary of couplings
    ## Output units: [eV]
    if cq["c4"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s"]:
            cp += (cq["c4"]*mp/m_q[q])*(m_q[q]-6*mBar)*Delta["p"][q] # factor of six is because I'm assuming we sum over all quarks
            cn += (cq["c4"]*mp/m_q[q])*(m_q[q]-6*mBar)*Delta["n"][q] # factor of six is because I'm assuming we sum over all quarks
        return [cp,cn]
    
def c5_N(cq):
    ## Define coupling, with the 1/Lambda^2 pulled out
    ## Input: dictionary of couplings
    ## Output units: unitless
    if cq["c5"]==0:
        return [0,0]
    else:
        cp = 3*cq["c5"]
        cn = 3*cq["c5"]
        return [cp,cn]

def c6_N(cq):
    ## Define coupling, with the 1/Lambda^2 pulled out
    ## Input: dictionary of couplings
    ## Output units: unitless
    if cq["c6"]==0:
        return [0,0]
    else:
        cp = 3*cq["c6"]
        cn = 3*cq["c6"]
        return [cp,cn]
    
def c7_N(cq):
    ## Define coupling, with the 1/Lambda^2 pulled out
    ## Input: dictionary of couplings
    ## Output units: unitless
    if cq["c7"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s","c","b","t"]:
            cp += cq["c7"]*Delta["p"][q]
            cn += cq["c7"]*Delta["n"][q]
        return [cp,cn]
    
def c8_N(cq):
    ## Define coupling, with the 1/Lambda^2 pulled out
    ## Input: dictionary of couplings
    ## Output units: unitless
    if cq["c8"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s","c","b","t"]:
            cp += cq["c8"]*Delta["p"][q]
            cn += cq["c8"]*Delta["n"][q]
        return [cp,cn]
    
def c9_N(cq):
    ## Define coupling, with the 1/Lambda^2 pulled out
    ## Input: dictionary of couplings
    ## Output units: unitless
    if cq["c9"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s","c","b","t"]:
            cp += cq["c9"]*delta["p"][q]
            cn += cq["c9"]*delta["n"][q]
        return [cp,cn]
    
def c10_N(cq):
    ## Define coupling, with the 1/Lambda^2 pulled out
    ## Input: dictionary of couplings
    ## Output units: unitless
    if cq["c10"]==0:
        return [0,0]
    else:
        cp = 0
        cn = 0
        for q in ["u","d","s","c","b","t"]:
            cp += cq["c10"]*delta["p"][q]
            cn += cq["c10"]*delta["n"][q]
        return [cp,cn]
    

######################################################################################################
## Map the relativistic couplings to a cross section (to give approx measure of interaction strength)
######################################################################################################
def sigma_from_EFT(cq,mX,Lambda):
    """
    [cq] = unitless set of active relativistic operators
    [mX] = [eV] DM mass
    [Lambda] = [eV] new physics scale

    Output units: [cm]^2 cross section of the nucleon 'vector'
    """
    c1 = c1_N(cq)
    c2 = c2_N(cq)
    c3 = c3_N(cq)
    c4 = c4_N(cq)
    c5 = c5_N(cq)
    c6 = c6_N(cq)
    c7 = c7_N(cq)
    c8 = c8_N(cq)
    c9 = c9_N(cq)
    c10 = c10_N(cq)
    c_sq = 0 # sum of c^2 values. Will have units of [eV]^-4
    ## sum for c1_NR
    c_sq+=(c1[0]/pow(Lambda,3)+c5[0]/pow(Lambda,2))**2
    c_sq+=(c1[1]/pow(Lambda,3)+c5[1]/pow(Lambda,2))**2
    ## sum for c4_NR
    c_sq+=((8*c9[0]-4*c8[0])/pow(Lambda,2))**2
    c_sq+=((8*c9[1]-4*c8[1])/pow(Lambda,2))**2
    ## sum for c6_NR
    c_sq+=((mp/mX)*c4[0]/pow(Lambda,3))**2
    c_sq+=((mp/mX)*c4[1]/pow(Lambda,3))**2
    ## sum for c7_NR
    c_sq+=(2*c7[0]/pow(Lambda,2))**2
    c_sq+=(2*c7[1]/pow(Lambda,2))**2
    ## sum for c8_NR
    c_sq+=(2*c6[0]/pow(Lambda,2))**2
    c_sq+=(2*c6[1]/pow(Lambda,2))**2
    ## sum for c9_NR
    c_sq+= (2*c6[0]/pow(Lambda,2)+2*(mp/mX)*c7[0]/pow(Lambda,2))**2
    c_sq+= (2*c6[1]/pow(Lambda,2)+2*(mp/mX)*c7[1]/pow(Lambda,2))**2
    ## sum for c10_NR
    c_sq+= (c3[0]/pow(Lambda,3)-2*(mp/mX)*c10[0]/pow(Lambda,2))**2
    c_sq+= (c3[1]/pow(Lambda,3)-2*(mp/mX)*c10[1]/pow(Lambda,2))**2
    ## sum for c11_NR
    c_sq+= (2*c10[0]/pow(Lambda,2)-2*(mp/mX)*c2[0]/pow(Lambda,3))**2
    c_sq+= (2*c10[1]/pow(Lambda,2)-2*(mp/mX)*c2[1]/pow(Lambda,3))**2
    ## compute cross section based on this sum
    sig = eV2_to_cm2*pow(mp*mX/(mp+mX),2)*c_sq/np.pi
    return sig


def sigma_p_from_EFT(cq,mX,Lambda):
    """
    [cq] = unitless set of active relativistic operators
    [mX] = [eV] DM mass
    [Lambda] = [eV] new physics scale

    Output units: [cm]^2 proton-DM cross section 
    """
    c1 = c1_N(cq)
    c2 = c2_N(cq)
    c3 = c3_N(cq)
    c4 = c4_N(cq)
    c5 = c5_N(cq)
    c6 = c6_N(cq)
    c7 = c7_N(cq)
    c8 = c8_N(cq)
    c9 = c9_N(cq)
    c10 = c10_N(cq)
    c_sq = 0 # sum of c^2 values. Will have units of [eV]^-4
    ## sum for c1_NR
    c_sq+=(c1[0]/pow(Lambda,3)+c5[0]/pow(Lambda,2))**2
    ## sum for c4_NR
    c_sq+=((8*c9[0]-4*c8[0])/pow(Lambda,2))**2
    ## sum for c6_NR
    c_sq+=((mp/mX)*c4[0]/pow(Lambda,3))**2
    ## sum for c7_NR
    c_sq+=(2*c7[0]/pow(Lambda,2))**2
    ## sum for c8_NR
    c_sq+=(2*c6[0]/pow(Lambda,2))**2
    ## sum for c9_NR
    c_sq+= (2*c6[0]/pow(Lambda,2)+2*(mp/mX)*c7[0]/pow(Lambda,2))**2
    ## sum for c10_NR
    c_sq+= (c3[0]/pow(Lambda,3)-2*(mp/mX)*c10[0]/pow(Lambda,2))**2
    ## sum for c11_NR
    c_sq+= (2*c10[0]/pow(Lambda,2)-2*(mp/mX)*c2[0]/pow(Lambda,3))**2
    ## compute cross section based on this sum
    sig = eV2_to_cm2*pow(mp*mX/(mp+mX),2)*c_sq/np.pi
    return sig

def sigma_n_from_EFT(cq,mX,Lambda):
    """
    [cq] = unitless set of active relativistic operators
    [mX] = [eV] DM mass
    [Lambda] = [eV] new physics scale

    Output units: [cm]^2 neutron-DM cross section 
    """
    c1 = c1_N(cq)
    c2 = c2_N(cq)
    c3 = c3_N(cq)
    c4 = c4_N(cq)
    c5 = c5_N(cq)
    c6 = c6_N(cq)
    c7 = c7_N(cq)
    c8 = c8_N(cq)
    c9 = c9_N(cq)
    c10 = c10_N(cq)
    c_sq = 0 # sum of c^2 values. Will have units of [eV]^-4
    ## sum for c1_NR
    c_sq+=(c1[1]/pow(Lambda,3)+c5[1]/pow(Lambda,2))**2
    ## sum for c4_NR
    c_sq+=((8*c9[1]-4*c8[1])/pow(Lambda,2))**2
    ## sum for c6_NR
    c_sq+=((mp/mX)*c4[1]/pow(Lambda,3))**2
    ## sum for c7_NR
    c_sq+=(2*c7[1]/pow(Lambda,2))**2
    ## sum for c8_NR
    c_sq+=(2*c6[1]/pow(Lambda,2))**2
    ## sum for c9_NR
    c_sq+= (2*c6[1]/pow(Lambda,2)+2*(mp/mX)*c7[1]/pow(Lambda,2))**2
    ## sum for c10_NR
    c_sq+= (c3[1]/pow(Lambda,3)-2*(mp/mX)*c10[1]/pow(Lambda,2))**2
    ## sum for c11_NR
    c_sq+= (2*c10[1]/pow(Lambda,2)-2*(mp/mX)*c2[1]/pow(Lambda,3))**2
    ## compute cross section based on this sum
    sig = eV2_to_cm2*pow(mp*mX/(mp+mX),2)*c_sq/np.pi
    return sig
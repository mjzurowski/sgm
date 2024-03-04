eV = 1 # base mass units
keV = 1E3 
MeV = 1E6 
GeV = 1E9

c = 299792458. # speed of light [km/s]
mp =  0.9314941*GeV # unified atomic mass
cpd_conversion = 86400*c*1E14 # for conversion of dsigdER to units of cpd/kg/keV
kg_to_eV = 1.8E-36 # conversion between eV and kg

Gf = 1.166364*1e-5 ##Fermi constant in natural unit GeV
sin2thetaW = 0.2229 ## Weak mixing angle from https://physics.nist.gov/cgi-bin/cuu/Value?sin2th
mproton=938.272*MeV
mn=939.565*MeV
mC12=12*mp
Delta_np = (mn-mproton)
M_nucleus=(mproton+mn)/2
me= 0.511*MeV
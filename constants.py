eV = 1 # base mass units
keV = 1E3 
GeV = 1E9

c = 299792458. # speed of light [m/s]
kms = c*1E-3 # get units of [km/s]
mp =  0.9314941*GeV # mass of a nucleon in eV.
per_day = 86400*c*1E2 # get from per cm to per day
cpd_conversion = per_day*keV*GeV # for conversion of dsigdER to units of cpd/kg/keV (given energy and rho units)
kg_to_eV = 1.8E-36 # conversion between eV and kg
eV2_to_cm2 = (197.236*1E-7)*(197.236*1E-7) # to get from units of inverse eV2 to cm2 (needed for switching between constructions of dsigdER)
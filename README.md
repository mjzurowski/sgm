## Signal Generation and Modelling

This code can be used to produce a DM event rate spectrum as a function of observed energy. A minimal example is shown in Demo.ipynb. The idea behind this is that a large number of different models/detectors/velocity distributions can be saved internally and combined in whatever way the user chooses. New versions of each can easily be included where desired, and used with existing versions of the others.
Further details about how the calculation is performed can be found in the pdf document attached (rate_calc)

The calculation is performed using a number of user defined objects:

1. Target (target.py)
A target is a nucleus able to undergo scattering. It should be defined with its mass number and nuclear form factors.
The current targets are available:
- Na (na.py)

2. DM model (dmmodel.py)
The user must define a DM interaction model (as a function of recoil energy) by creating a new "DMModel" object which requires definition of a minimum velocity (in units of km/s) and a differential cross section (in units of cm2/eV). These will both (probably) be a function of target objects which will be read in as an argument for the functions.
The current models are available (see wimp.py):
- Standard SI WIMP (computed using Helm form factor)
- NREFT WIMP defined with operators O1, O4, O5, O6, O7.

3. Detector (detector.py)
A detector object has associated with it some number of targets. The easiest way of doing this is to define the targets separately, then import them into whereever you define your detector. The user needs to define these "Nuclei" as a list, along with the associated transformation between recoil energy and observed energy "ER_R" (e.g., quenching factor, ionisation yield). The energy region of interest (ROI), max energy, resolution functions (DeltaE, Res) and efficiency should all also be defined as a function of observed energy.
The current detectors are available:
- DAMA (nai.py)

4. Velocity distribution (veldists.py)
Velocity distributions should be saved as .dat files in the velocity_distributions folder. They should be given as the total velocity integral computed as a function of minimum velocity. These can then be read in with an associated DM density.
The current distributions are available:
- SHM
- SHM++
- Shards
- Streams
- Lawrence et al. (8 different distributions)

### Use
A basic example of how this code can be used is shown in Demo.ipynb with detectors/models/velocity distributions that already exist. There are a few required steps:
1. Initialise your detector, model, and velocity distribution.
2. Choose a particular DM mass and cross section, and use these to define a DM interaction rate as a function of target (T) and recoil energy (E) for the Model in use with Model.dRdER.
3. Compute this rate as a function of observed energy for the detector by passing it as an argument in Det.dRdR.
## Signal Generation and Modelling

This code can be used to produce a DM event rate spectrum as a function of observed energy. A minimal example is shown in Demo.ipynb. The idea behind this is that a large number of different models/detectors/velocity distributions can be saved internally and combined in whatever way the user chooses. New versions of each can easily be included where desired, and used with existing versions of the others.
Further details about how the calculation is performed can be found in the pdf document attached (rate_calc).

The calculation is performed using a number of user defined objects:

1. Target (target.py)
A target is a nucleus able to undergo scattering. It should be defined with its mass number, nuclear form factors, electron transition probabilities, photoelectric absorption. Nuclear form factors are assumed to be defined in the same form as in Fitzpatrick et al. (https://arxiv.org/abs/1203.3542)
The current targets are available:
- Na (na.py)
- I (i.py)

2. DM model (dmmodel.py)
The user must define a DM interaction model (as a function of recoil energy) by creating a new "DMModel" object which requires definition of a minimum velocity (in units of km/s) and a differential cross section (in units of cm2/eV). These will both (probably) be a function of target objects which will be read in as an argument for the functions.
The current models are available (see models/):
- Standard SI WIMP (computed using Helm form factor)
- NREFT WIMP defined with operators O1, O4, O5, O6, O7 (see appendix A.3 of https://arxiv.org/abs/1203.3542).
- SI inelastic WIMP (an example of a model that has requires additional variables)

3. Detector (detector.py)
A detector object has associated with it some number of targets. The easiest way of doing this is to define the targets separately, then import them into whereever you define your detector. The user needs to define these "Nuclei" as a list, along with the associated transformation between recoil energy and observed energy "ER_R" (e.g., quenching factor, ionisation yield). The energy region of interest (ROI), max energy, resolution functions (DeltaE, Res) and efficiency should all also be defined as a function of observed energy.
If this is being used with an ER interaction model (i.e., backgrounds), dRdE can be called with the flag NR=False and this will prevent the NR --> observed signal transformations from being applied.
The current detectors are available (see detectors/):
- DAMA (dama.py)
- COSINE (cosine.py)

4. Velocity distribution (veldists.py)
Velocity distributions should be saved as .dat files in the velocity_distributions folder. They should be given as the total velocity integral computed as a function of minimum velocity. These can then be read in with an associated DM density. See the rate_calc document for a more detailed explanation of the format of these distributions, and references to the distributions defined below.
The current distributions are available:
- SHM
- SHM++ (https://arxiv.org/abs/1810.11468)
- Shards (https://arxiv.org/abs/1909.04684)
- Streams (https://arxiv.org/abs/1807.09004)
- Lawrence et al. (8 different distributions) (https://arxiv.org/abs/2207.07644)

### Use
A basic example of how this code can be used is shown in Demo.ipynb with detectors/models/velocity distributions that already exist. There are a few required steps:
1. Initialise your detector, model, and velocity distribution.
2. Call Det.dRdE, which will calculate the rate observed in the detector. This takes as variables the energy, the model (Model.dRdER) and the DM parameters the model requires, passed as kwargs (e.g., mass, cross section, velocity distribution etc)

### To do
- ~Test adding detectors with more complex quenching factors~ Added using simple interpolation rather than fitting to a proper model
- Optimise integral in detector rate calc
- ~Deal with models that have g and h velocity integrals better~ Distribution objects now have both h and g integrals that can be called separately when defining model rates
- ~Improve final "output" - probably want it to be a function of mX and sig~ Added kwargs to the detector functions to allow more generic rate forms/variables to be passed directly to Det.dRdE
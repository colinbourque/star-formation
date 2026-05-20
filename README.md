# star-formation
Files, scripts, and assorted data used in the radiative transfer star formation modeling project. For all functions pertaining to the generation of a model, and mathematical analysis, the only dependencies are NumPy as Astropy. Plotting, if any scripts for that end up being posted here, is all built around Matplotlib.

/make_files/ contains files needed for generating runs of radmc3d, including the underlying physical dust density profiles and simulation grids. .inp type files in the main directory here are inhereted from the 2d run of the code. These are then handled by the .py files and turned into radmc3d compatible input files, which are stored in /make_files/aux/  

/analysis/ contains script(s) for interpreting the outputs of radmc3d

/model_parameters/ contains files describing the inputs to the simulations

/scripts/ is currently only in use for passing .sh files back and forth between my local machine and the computational tower, since it is easier for me to generate them locally.

## Running the code

#### Pre-requisites

Building the model data requires NumPy (>=1.26.4) and Astropy (>=6.0.1), and running the radiative transfer code requires an installation of radmc3d. Tentatively, this codebase is written to only depend on .py scripts; although if I do include any notebooks here, running them will require Jupyter.

#### Directory structure

This repository contains all of the scripts needed to generate input files, control radmc3d, and analyze outputs. For the actual handling of these input and output files, I operate in a separate directory which I will call /simdir/. This directory contains a copy of the model parameters tables, as well as folders for each model's respective input/output/analysis files. A visualization of this file structure is attached here:

simdir/
├── model_parameters/
│   └── model_parameters_1.tbl
└── model01/
    ├── AUX
    ├── INPUT
    ├── RUNDIR
    ├── OUTPUT
    └── RESULTS
star-formation/
└── make_files/
    ├── problem_setup.py
    ├── radmc_utils.py
    ├── dustopac_1.inp
    ├── external_meanint.inp
    ├── frequency.inp
    └── grid.plt34.mod
└── analysis/
└── scripts/
└── model_parameters/

Most scripts rely on directory pointings which are coded relative to /simdir/, via a variable $\texttt{simdir}$ which should be re-written to point towards whatever directory this is on the machine where the code is being run. The /star-formation/ directory is set up to match this repository, and its own pointings likely do not require modification.

In order to begin a full model run, it is necessary to create a corresponding directory and subdirectories within /simdir/. While /simdir/ can be anything, as of right now this code is written to explicitly look for an internal structure which matches that shown above.

#### Making the input and auxilliary files

Generation of the model input files and radmc3d control files is done via make_files/problem_setup.py. After making sure to set $\texttt{simdir}$ correctly, run this script

>python problem_setup.py

You will then be prompted for the model number you'd like to generate. As of right now, this has protections coded in to ensure an integer is provided (and therefore avoid a litany of errors) however this could be removed to allow for more granular model naming (eg variations could be named 1a, 1b, etc). Upon providing the model number, this will generate a set of radmc3d control files to /simdir/modelNN/AUX/ and the grid, density, and stellar source input files to /simdir/modelNN/INPUT/. 

#### Running radmc3d

Unfortunately, radmc3d is quite rigid in the naming conventions it expects for input files. As such, a secondary temporary directory, /RUNDIR/ is established for actually running the radiative transfer code. 

For every timestep, the following must occur:
0. The auxillairy files are copied to RUNDIR
1. stars.NNN.inp, amr_grid.NNN.inp, and dust_density.NNN.inp are copied to RUNDIR without the model number .NNN label.
2. Thermal radiative transfer is computed, and the temperatures file is copied to /OUTPUT/
3. The ISRF photons (external_source.inp) are disabled (renamed to something that radmc3d won't find)
4. Spectra for all inclinations are generated and copied to /OUTPUT/ with modelnumber and inclination labels
5. The ISRF photons are turned back on (a new external_source.inp is copied over from AUX)
6. Step up in model number and repeat beginning with step 1.

Inside of /scripts/ I have attached a Python script (runscript.py) which will generate a file (runmodNN.sh) in /RUNDIR/ that can be used to run this process in bulk. Again, presuming /simdir/ is set correctly, this script can be run as

>python runscript.py

with the model number provided subsequently. 








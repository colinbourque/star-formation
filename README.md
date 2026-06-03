# star-formation
Files, scripts, and assorted data used in the radiative transfer star formation modeling project. For all functions pertaining to the generation of a model, and mathematical analysis, the only dependencies are `NumPy` and `Astropy`. Plotting, if any scripts for that end up being posted here, is all built around `Matplotlib`.

`/make_files/` contains files needed for generating runs of `radmc3d`, including the underlying physical dust density profiles and simulation grids. `.inp` type files in the main directory here are inhereted from the 2d run of the code. These are then handled by the `.py` files and turned into `radmc3d` compatible input files, which are stored in each respective model's `/simdir/AUX/` directory.  

`/analysis/` contains script(s) for interpreting the outputs of `radmc3d`

`/model_parameters/` contains files describing the inputs to the simulations

`/rad-transfer/` contains scripts used for simplifying the process of actually running radmc3d

`/scripts/` is currently only in use for passing .sh files back and forth between my local machine and the computational tower, since it is easier for me to generate them locally.

`/INDICATORS/` contains compiled results from each model run, both bolometric luminosity and temperature tables, and figures comparing the 2D results to the 3D results.

## More detailed documentation

The remainder of the `readme` file contains a cursory explanation of how the code is meant to be run. For a more detailed explanation of the physical workings of the code please see the attached document `star-formation-docs-v2.pdf`.

## Running the code

#### Pre-requisites

Building the model data requires NumPy (>=1.26.4) and Astropy (>=6.0.1), and running the radiative transfer code requires an installation of radmc3d. Tentatively, this codebase is written to only depend on .py scripts; although if I do include any notebooks here, running them will require Jupyter.

#### Directory structure

This repository contains all of the scripts needed to generate input files, control radmc3d, and analyze outputs. For the actual handling of these input and output files, I operate in a separate directory which I will call /simdir/. This directory contains a copy of the model parameters tables, as well as folders for each model's respective input/output/analysis files. A visualization of this file structure is attached here: 

<img width="311" height="400" alt="Screenshot 2026-05-20 at 1 47 13 PM" src="https://github.com/user-attachments/assets/bb9d5f28-de62-4ce5-9767-376230785c4d" />


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

Inside of /rad-transfer/ I have attached a Python script (make_runscript.py) which will generate a file (runmodNN.sh) in /RUNDIR/ that can be used to run this process in bulk. Again, presuming /simdir/ is set correctly, this script can be run as

>python runscript.py

with the model number provided subsequently. This will generate a script runmodNN.sh inside of /RUNDIR/ which can then be run to make radmc3d compute all products for a given model.

As an alternative method for generating more complex configurations of runscripts, a notebook `/scripts/runscript.ipynb` is also included. Here, an argument for number of partitions can be set (in addition to the standard `simdir` and `nthreads`). This will create multiple `runmodNN-P.sh` scripts containing commands to run \[number of timesteps\]//p timesteps, labeled "A" through \["A", "B", "C", ...]\[P\]. <ins>Each partition</ins> will use `nthreads` cores, meaning to run all partitions simtaneously, at least nthreads*P cores should be available.

#### Getting bolometric temperatures and luminosities

Bolometric temperatures and luminosities are calculated by integrating the output spectrum from `radmc3d` (using function `trapezoid()`. Functions `bol_luminosity(nu, fnu)` and `bol_temperature(nu, fnu)` are included in `radmc_utils.py`, and return values in units of $L_\odot$ at 1 PC, and K, respectively.

For convenience, a script `calc_lbol_tbol.py` is included in `/analysis/`. Pre-requisites for running this script are to ensure that `simdir` is set correctly, and that this script has access to `radmc_utils.py`. The easiest way to achieve this latter pre-requisite is by just copying one or the other into the same folder (the difficult way would be to turn `radmc_utils.py` into a package and install it to your working environment). Presuming this is all set, then running 

>python get_lbol_tbol.py

will prompt for a model number, and, when provided one, will create two files `lbols_modelNN.tab` and `tbols_modelNN.tab` in directory `./indicators/`. Both of these tables are structures to contain the time in Myr, $\Delta$ time in yr, and then eight columns of indicators progressing from $5\degree$ indlication to $85\degree$. These labels are written into a header row, as well.

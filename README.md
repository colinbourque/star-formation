# star-formation
Files, scripts, and assorted data used in the radiative transfer star formation modeling project. For all functions pertaining to the generation of a model, and mathematical analysis, the only dependency is NumPy. Plotting, if any scripts for that end up being posted here, is all built around Matplotlib.

/make_files/ contains files needed for generating runs of radmc3d, including the underlying physical dust density profiles and simulation grids. .inp type files in the main directory here are inhereted from the 2d run of the code. These are then handled by the .py files and turned into radmc3d compatible input files, which are stored in /make_files/aux/  

/model_parameters/ contains files describing the inputs to the simulations

/scripts/ is currently only in use for passing .sh files back and forth between my local machine and the computational tower, since it is easier for me to generate them locally.

import numpy as np

print(f'Make sure variable simdir is set correctly.')
modelnumber = input('Generate run script for model number...')

mdir   = f'/Users/cbourque/astrophysics/star-formation/SIMULATIONS/model{modelnumber:02n}'
auxdir = f'{mdir}/AUX'
indir  = f'{mdir}/INPUT'
outdir = f'{mdir}/OUTPUT'
rundir = f'{mdir}/RUNDIR'

mpars = np.loadtxt(f'/Users/colinbourque/astrophysics/star-formation/SIMULATIONS/model_parameters/model_parameters_{modelnumber}.tbl', skiprows=1)
mpars = mpars[~np.isnan(mpars).any(axis=1)].T ## remove entries with nan values

nthreads = 4

with open(f'{rundir}/runmod{modelnumber}.sh','w') as f:
    ## copy the auxilliary files into the working 'run' directory
    f.write(f'export STARTTIME=$(date) \n')
    f.write(f'cp {mdir}/AUX/aperture_info.inp {rundir}/aperture_info.inp \n')
    f.write(f'cp {mdir}/AUX/dustopac.inp {rundir}/dustopac.inp \n')
    f.write(f'cp {mdir}/AUX/radmc3d.inp {rundir}/radmc3d.inp \n')
    f.write(f'cp {mdir}/AUX/wavelength_micron.inp {rundir}/wavelength_micron.inp \n')
    f.write(f'cp {mdir}/AUX/dustkappa_OH5.inp {rundir}/dustkappa_OH5.inp \n')

    for i in range(1,len(mpars[0])): ## call this to run everything at once
        f.write(f'echo --------------------------------------------------------- \n')
        f.write(f'echo BEGIN MODEL STEP {i} \n')
        f.write(f'echo --------------------------------------------------------- \n')

        ## copy each timestep's files into the run directory
        f.write(f'cp {mdir}/INPUT/stars.{i:03n}.inp {rundir}/stars.inp \n')
        f.write(f'cp {mdir}/INPUT/amr_grid.{i:03n}.inp {rundir}/amr_grid.inp \n')
        f.write(f'cp {mdir}/INPUT/dust_density.{i:03n}.inp {rundir}/dust_density.inp \n')

        ## run the thermal monte carlo
        f.write(f'cp {mdir}/AUX/external_source.inp {rundir}/external_source.inp \n')
        f.write(f'radmc3d mctherm setthreads {nthreads} countwrite 1000000 \n')
        f.write(f'cp {rundir}/dust_temperature.dat {mdir}/OUTPUT/dust_temperature.{i:03n}.dat \n')
        f.write(f'mv {rundir}/external_source.inp {rundir}/no_external_source.inp \n')

        ## make SEDs
        for incl in [5,15,25,35,45,55,65,75,85]:
            f.write(f'radmc3d sed incl {incl} useapert dpc 140 setthreads {nthreads}\n')
            f.write(f'mv {rundir}/spectrum.out {mdir}/OUTPUT/spectrum.{i:03n}.incl{incl:02n}.out \n')

    f.write(f'echo STARTED AT $STARTTIME \n')
    f.write(f'echo ENDED AT $(date) \n')

export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 391 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.391.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.391.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.391.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.391.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.391.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 392 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.392.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.392.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.392.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.392.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.392.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 393 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.393.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.393.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.393.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.393.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.393.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 394 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.394.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.394.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.394.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.394.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.394.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 395 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.395.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.395.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.395.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.395.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.395.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 396 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.396.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.396.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.396.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.396.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.396.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 397 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.397.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.397.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.397.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.397.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.397.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 398 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.398.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.398.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.398.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.398.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.398.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 399 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.399.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.399.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.399.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.399.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.399.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 400 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.400.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.400.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.400.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.400.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.400.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 401 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.401.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.401.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.401.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.401.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.401.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 402 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.402.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.402.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.402.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.402.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.402.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 403 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.403.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.403.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.403.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.403.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.403.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 404 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.404.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.404.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.404.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.404.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.404.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 405 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.405.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.405.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.405.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.405.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.405.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 406 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.406.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.406.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.406.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.406.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.406.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 407 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.407.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.407.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.407.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.407.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.407.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 408 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.408.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.408.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.408.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.408.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.408.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 409 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.409.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.409.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.409.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.409.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.409.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 410 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.410.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.410.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.410.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.410.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.410.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 411 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.411.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.411.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.411.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.411.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.411.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 412 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.412.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.412.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.412.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.412.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.412.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 413 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.413.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.413.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.413.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.413.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.413.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 414 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.414.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.414.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.414.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.414.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.414.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 415 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.415.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.415.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.415.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.415.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.415.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 416 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.416.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.416.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.416.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.416.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.416.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 417 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.417.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.417.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.417.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.417.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.417.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 418 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.418.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.418.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.418.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.418.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.418.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 419 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.419.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.419.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.419.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.419.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.419.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 420 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.420.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.420.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.420.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.420.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.420.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 421 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.421.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.421.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.421.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.421.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.421.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 422 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.422.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.422.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.422.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.422.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.422.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 423 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.423.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.423.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.423.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.423.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.423.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 424 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.424.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.424.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.424.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.424.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.424.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 425 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.425.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.425.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.425.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.425.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.425.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 426 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.426.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.426.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.426.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.426.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.426.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 427 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.427.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.427.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.427.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.427.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.427.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 428 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.428.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.428.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.428.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.428.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.428.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 429 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.429.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.429.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.429.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.429.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.429.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 430 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.430.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.430.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.430.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.430.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.430.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 431 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.431.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.431.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.431.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.431.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.431.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 432 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.432.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.432.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.432.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.432.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.432.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 433 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.433.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.433.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.433.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.433.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.433.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 434 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.434.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.434.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.434.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.434.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.434.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 435 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.435.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.435.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.435.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.435.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.435.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 436 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.436.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.436.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.436.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.436.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.436.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 437 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.437.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.437.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.437.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.437.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.437.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 438 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.438.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.438.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.438.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.438.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.438.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 439 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.439.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.439.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.439.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.439.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.439.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 440 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.440.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.440.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.440.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.440.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.440.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 441 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.441.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.441.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.441.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.441.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.441.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 442 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.442.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.442.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.442.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.442.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.442.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 443 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.443.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.443.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.443.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.443.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.443.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 444 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.444.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.444.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.444.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.444.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.444.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 445 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.445.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.445.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.445.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.445.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.445.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 446 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.446.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.446.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.446.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.446.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.446.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 447 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.447.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.447.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.447.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.447.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.447.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 448 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.448.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.448.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.448.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.448.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.448.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 449 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.449.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.449.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.449.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.449.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.449.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 450 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.450.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.450.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.450.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.450.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.450.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 451 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.451.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.451.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.451.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.451.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.451.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 452 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.452.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.452.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.452.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.452.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.452.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 453 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.453.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.453.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.453.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.453.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.453.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 454 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.454.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.454.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.454.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.454.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.454.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 455 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.455.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.455.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.455.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.455.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.455.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 456 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.456.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.456.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.456.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.456.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.456.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 457 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.457.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.457.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.457.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.457.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.457.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 458 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.458.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.458.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.458.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.458.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.458.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 459 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.459.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.459.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.459.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.459.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.459.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 460 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.460.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.460.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.460.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.460.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.460.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 461 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.461.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.461.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.461.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.461.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.461.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 462 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.462.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.462.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.462.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.462.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.462.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 463 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.463.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.463.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.463.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.463.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.463.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 464 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.464.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.464.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.464.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.464.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.464.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 465 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.465.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.465.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.465.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.465.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.465.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 466 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.466.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.466.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.466.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.466.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.466.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 467 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.467.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.467.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.467.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.467.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.467.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 468 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.468.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.468.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.468.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.468.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.468.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 469 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.469.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.469.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.469.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.469.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.469.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 470 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.470.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.470.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.470.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.470.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.470.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 471 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.471.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.471.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.471.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.471.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.471.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 472 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.472.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.472.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.472.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.472.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.472.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 473 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.473.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.473.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.473.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.473.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.473.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 474 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.474.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.474.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.474.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.474.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.474.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 475 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.475.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.475.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.475.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.475.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.475.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 476 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.476.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.476.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.476.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.476.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.476.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 477 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.477.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.477.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.477.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.477.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.477.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 478 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.478.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.478.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.478.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.478.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.478.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 479 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.479.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.479.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.479.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.479.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.479.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 480 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.480.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.480.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.480.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.480.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.480.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 481 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.481.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.481.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.481.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.481.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.481.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 482 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.482.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.482.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.482.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.482.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.482.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 483 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.483.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.483.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.483.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.483.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.483.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 484 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.484.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.484.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.484.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.484.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.484.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 485 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.485.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.485.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.485.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.485.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.485.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 486 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.486.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.486.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.486.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.486.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.486.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 487 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.487.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.487.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.487.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.487.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.487.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 488 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.488.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.488.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.488.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.488.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.488.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 489 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.489.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.489.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.489.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.489.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.489.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 490 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.490.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.490.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.490.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.490.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.490.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 491 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.491.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.491.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.491.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.491.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.491.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 492 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.492.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.492.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.492.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.492.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.492.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 493 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.493.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.493.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.493.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.493.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.493.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 494 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.494.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.494.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.494.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.494.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.494.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 495 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.495.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.495.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.495.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.495.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.495.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 496 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.496.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.496.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.496.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.496.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.496.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 497 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.497.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.497.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.497.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.497.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.497.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 498 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.498.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.498.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.498.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.498.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.498.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 499 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.499.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.499.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.499.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.499.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.499.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 500 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.500.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.500.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.500.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.500.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.500.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 501 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.501.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.501.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.501.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.501.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.501.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 502 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.502.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.502.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.502.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.502.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.502.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 503 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.503.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.503.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.503.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.503.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.503.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 504 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.504.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.504.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.504.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.504.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.504.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 505 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.505.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.505.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.505.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.505.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.505.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 506 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.506.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.506.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.506.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.506.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.506.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 507 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.507.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.507.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.507.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.507.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.507.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 508 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.508.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.508.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.508.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.508.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.508.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 509 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.509.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.509.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.509.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.509.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.509.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 510 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.510.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.510.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.510.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.510.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.510.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 511 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.511.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.511.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.511.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.511.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.511.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 512 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.512.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.512.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.512.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.512.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.512.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 513 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.513.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.513.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.513.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.513.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.513.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 514 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.514.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.514.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.514.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.514.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.514.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 515 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.515.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.515.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.515.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.515.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.515.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 516 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.516.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.516.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.516.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.516.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.516.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 517 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.517.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.517.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.517.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.517.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.517.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 518 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.518.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.518.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.518.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.518.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.518.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 519 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.519.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.519.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.519.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.519.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.519.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 520 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.520.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.520.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.520.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.520.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.520.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 521 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.521.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.521.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.521.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.521.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.521.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 522 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.522.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.522.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.522.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.522.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.522.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 523 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.523.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.523.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.523.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.523.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.523.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 524 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.524.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.524.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.524.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.524.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.524.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 525 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.525.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.525.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.525.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.525.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.525.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 526 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.526.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.526.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.526.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.526.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.526.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 527 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.527.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.527.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.527.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.527.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.527.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 528 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.528.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.528.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.528.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.528.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.528.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 529 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.529.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.529.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.529.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.529.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.529.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 530 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.530.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.530.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.530.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.530.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.530.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 531 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.531.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.531.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.531.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.531.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.531.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 532 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.532.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.532.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.532.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.532.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.532.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 533 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.533.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.533.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.533.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.533.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.533.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 534 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.534.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.534.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.534.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.534.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.534.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 535 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.535.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.535.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.535.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.535.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.535.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 536 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.536.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.536.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.536.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.536.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.536.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 537 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.537.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.537.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.537.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.537.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.537.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 538 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.538.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.538.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.538.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.538.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.538.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 539 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.539.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.539.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.539.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.539.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.539.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 540 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.540.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.540.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.540.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.540.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.540.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 541 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.541.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.541.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.541.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.541.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.541.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 542 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.542.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.542.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.542.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.542.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.542.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 543 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.543.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.543.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.543.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.543.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.543.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 544 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.544.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.544.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.544.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.544.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.544.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 545 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.545.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.545.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.545.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.545.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.545.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 546 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.546.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.546.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.546.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.546.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.546.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 547 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.547.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.547.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.547.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.547.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.547.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 548 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.548.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.548.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.548.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.548.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.548.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 549 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.549.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.549.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.549.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.549.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.549.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 550 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.550.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.550.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.550.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.550.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.550.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 551 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.551.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.551.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.551.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.551.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.551.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 552 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.552.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.552.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.552.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.552.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.552.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 553 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.553.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.553.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.553.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.553.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.553.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 554 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.554.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.554.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.554.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.554.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.554.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 555 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.555.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.555.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.555.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.555.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.555.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 556 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.556.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.556.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.556.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.556.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.556.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 557 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.557.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.557.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.557.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.557.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.557.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 558 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.558.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.558.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.558.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.558.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.558.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 559 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.559.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.559.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.559.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.559.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.559.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 560 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.560.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.560.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.560.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.560.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.560.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 561 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.561.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.561.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.561.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.561.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.561.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 562 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.562.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.562.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.562.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.562.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.562.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 563 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.563.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.563.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.563.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.563.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.563.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 564 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.564.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.564.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.564.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.564.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.564.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 565 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.565.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.565.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.565.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.565.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.565.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 566 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.566.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.566.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.566.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.566.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.566.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 567 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.567.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.567.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.567.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.567.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.567.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 568 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.568.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.568.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.568.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.568.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.568.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 569 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.569.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.569.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.569.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.569.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.569.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 570 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.570.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.570.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.570.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.570.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.570.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 571 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.571.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.571.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.571.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.571.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.571.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 572 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.572.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.572.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.572.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.572.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.572.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 573 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.573.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.573.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.573.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.573.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.573.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 574 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.574.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.574.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.574.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.574.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.574.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 575 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.575.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.575.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.575.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.575.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.575.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 576 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.576.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.576.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.576.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.576.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.576.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 577 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.577.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.577.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.577.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.577.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.577.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 578 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.578.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.578.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.578.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.578.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.578.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 579 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.579.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.579.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.579.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.579.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.579.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 580 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.580.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.580.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.580.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.580.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.580.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

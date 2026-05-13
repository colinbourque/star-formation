export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 1 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.001.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.001.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.001.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.001.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.001.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 2 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.002.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.002.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.002.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.002.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.002.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 3 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.003.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.003.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.003.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.003.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.003.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 4 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.004.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.004.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.004.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.004.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.004.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 5 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.005.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.005.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.005.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.005.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.005.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 6 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.006.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.006.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.006.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.006.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.006.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 7 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.007.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.007.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.007.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.007.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.007.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 8 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.008.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.008.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.008.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.008.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.008.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 9 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.009.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.009.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.009.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.009.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.009.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 10 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.010.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.010.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.010.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.010.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.010.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 11 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.011.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.011.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.011.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.011.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.011.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 12 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.012.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.012.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.012.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.012.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.012.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 13 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.013.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.013.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.013.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.013.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.013.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 14 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.014.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.014.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.014.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.014.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.014.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 15 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.015.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.015.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.015.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.015.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.015.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 16 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.016.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.016.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.016.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.016.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.016.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 17 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.017.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.017.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.017.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.017.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.017.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 18 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.018.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.018.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.018.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.018.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.018.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 19 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.019.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.019.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.019.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.019.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.019.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 20 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.020.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.020.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.020.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.020.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.020.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 21 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.021.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.021.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.021.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.021.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.021.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 22 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.022.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.022.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.022.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.022.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.022.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 23 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.023.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.023.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.023.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.023.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.023.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 24 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.024.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.024.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.024.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.024.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.024.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 25 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.025.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.025.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.025.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.025.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.025.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 26 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.026.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.026.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.026.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.026.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.026.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 27 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.027.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.027.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.027.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.027.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.027.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 28 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.028.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.028.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.028.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.028.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.028.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 29 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.029.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.029.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.029.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.029.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.029.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 30 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.030.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.030.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.030.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.030.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.030.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 31 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.031.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.031.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.031.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.031.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.031.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 32 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.032.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.032.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.032.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.032.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.032.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 33 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.033.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.033.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.033.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.033.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.033.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 34 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.034.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.034.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.034.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.034.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.034.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 35 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.035.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.035.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.035.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.035.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.035.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 36 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.036.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.036.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.036.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.036.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.036.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 37 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.037.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.037.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.037.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.037.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.037.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 38 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.038.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.038.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.038.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.038.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.038.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 39 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.039.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.039.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.039.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.039.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.039.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 40 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.040.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.040.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.040.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.040.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.040.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 41 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.041.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.041.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.041.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.041.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.041.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 42 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.042.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.042.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.042.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.042.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.042.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 43 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.043.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.043.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.043.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.043.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.043.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 44 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.044.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.044.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.044.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.044.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.044.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 45 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.045.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.045.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.045.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.045.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.045.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 46 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.046.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.046.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.046.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.046.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.046.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 47 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.047.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.047.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.047.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.047.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.047.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 48 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.048.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.048.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.048.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.048.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.048.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 49 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.049.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.049.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.049.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.049.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.049.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 50 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.050.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.050.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.050.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.050.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.050.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 51 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.051.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.051.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.051.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.051.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.051.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 52 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.052.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.052.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.052.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.052.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.052.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 53 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.053.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.053.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.053.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.053.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.053.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 54 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.054.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.054.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.054.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.054.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.054.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 55 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.055.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.055.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.055.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.055.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.055.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 56 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.056.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.056.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.056.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.056.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.056.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 57 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.057.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.057.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.057.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.057.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.057.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 58 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.058.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.058.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.058.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.058.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.058.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 59 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.059.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.059.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.059.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.059.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.059.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 60 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.060.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.060.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.060.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.060.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.060.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 61 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.061.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.061.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.061.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.061.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.061.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 62 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.062.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.062.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.062.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.062.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.062.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 63 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.063.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.063.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.063.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.063.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.063.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 64 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.064.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.064.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.064.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.064.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.064.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 65 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.065.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.065.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.065.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.065.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.065.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 66 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.066.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.066.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.066.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.066.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.066.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 67 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.067.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.067.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.067.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.067.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.067.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 68 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.068.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.068.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.068.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.068.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.068.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 69 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.069.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.069.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.069.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.069.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.069.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 70 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.070.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.070.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.070.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.070.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.070.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 71 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.071.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.071.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.071.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.071.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.071.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 72 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.072.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.072.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.072.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.072.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.072.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 73 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.073.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.073.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.073.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.073.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.073.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 74 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.074.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.074.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.074.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.074.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.074.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 75 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.075.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.075.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.075.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.075.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.075.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 76 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.076.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.076.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.076.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.076.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.076.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 77 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.077.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.077.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.077.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.077.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.077.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 78 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.078.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.078.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.078.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.078.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.078.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 79 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.079.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.079.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.079.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.079.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.079.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 80 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.080.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.080.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.080.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.080.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.080.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 81 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.081.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.081.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.081.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.081.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.081.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 82 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.082.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.082.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.082.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.082.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.082.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 83 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.083.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.083.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.083.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.083.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.083.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 84 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.084.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.084.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.084.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.084.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.084.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 85 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.085.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.085.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.085.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.085.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.085.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 86 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.086.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.086.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.086.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.086.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.086.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 87 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.087.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.087.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.087.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.087.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.087.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 88 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.088.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.088.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.088.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.088.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.088.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 89 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.089.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.089.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.089.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.089.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.089.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 90 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.090.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.090.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.090.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.090.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.090.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 91 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.091.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.091.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.091.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.091.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.091.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 92 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.092.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.092.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.092.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.092.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.092.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 93 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.093.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.093.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.093.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.093.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.093.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 94 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.094.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.094.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.094.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.094.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.094.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 95 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.095.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.095.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.095.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.095.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.095.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 96 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.096.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.096.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.096.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.096.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.096.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 97 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.097.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.097.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.097.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.097.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.097.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 98 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.098.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.098.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.098.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.098.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.098.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 99 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.099.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.099.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.099.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.099.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.099.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 100 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.100.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.100.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.100.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.100.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.100.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 101 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.101.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.101.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.101.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.101.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.101.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 102 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.102.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.102.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.102.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.102.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.102.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 103 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.103.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.103.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.103.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.103.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.103.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 104 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.104.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.104.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.104.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.104.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.104.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 105 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.105.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.105.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.105.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.105.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.105.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 106 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.106.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.106.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.106.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.106.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.106.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 107 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.107.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.107.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.107.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.107.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.107.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 108 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.108.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.108.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.108.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.108.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.108.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 109 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.109.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.109.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.109.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.109.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.109.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 110 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.110.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.110.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.110.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.110.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.110.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 111 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.111.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.111.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.111.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.111.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.111.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 112 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.112.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.112.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.112.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.112.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.112.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 113 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.113.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.113.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.113.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.113.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.113.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 114 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.114.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.114.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.114.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.114.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.114.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 115 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.115.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.115.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.115.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.115.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.115.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 116 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.116.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.116.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.116.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.116.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.116.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 117 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.117.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.117.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.117.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.117.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.117.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 118 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.118.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.118.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.118.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.118.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.118.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 119 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.119.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.119.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.119.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.119.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.119.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 120 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.120.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.120.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.120.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.120.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.120.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 121 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.121.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.121.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.121.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.121.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.121.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 122 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.122.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.122.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.122.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.122.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.122.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 123 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.123.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.123.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.123.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.123.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.123.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 124 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.124.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.124.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.124.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.124.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.124.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 125 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.125.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.125.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.125.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.125.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.125.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 126 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.126.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.126.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.126.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.126.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.126.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 127 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.127.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.127.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.127.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.127.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.127.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 128 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.128.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.128.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.128.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.128.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.128.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 129 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.129.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.129.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.129.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.129.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.129.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 130 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.130.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.130.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 131 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.131.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.131.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 132 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.132.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.132.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 133 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.133.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.133.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 134 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.134.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.134.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 135 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.135.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.135.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 136 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.136.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.136.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 137 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.137.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.137.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 138 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.138.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.138.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 139 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.139.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.139.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 140 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.140.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.140.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 141 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.141.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.141.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 142 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.142.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.142.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 143 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.143.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.143.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 144 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/stars.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/amr_grid.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/INPUT/dust_density.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/dust_temperature.144.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-A/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/OUTPUT/spectrum.144.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

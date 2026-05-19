export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model20/RUNDIR-D/dustkappa_OH5.inp 
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

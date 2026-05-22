export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 130 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.130.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.130.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 131 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.131.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.131.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 132 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.132.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.132.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 133 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.133.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.133.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 134 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.134.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.134.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 135 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.135.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.135.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 136 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.136.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.136.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 137 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.137.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.137.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 138 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.138.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.138.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 139 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.139.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.139.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 140 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.140.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.140.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 141 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.141.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.141.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 142 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.142.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.142.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 143 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.143.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.143.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 144 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.144.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.144.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 145 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.145.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.145.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.145.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.145.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.145.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 146 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.146.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.146.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.146.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.146.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.146.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 147 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.147.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.147.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.147.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.147.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.147.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 148 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.148.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.148.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.148.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.148.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.148.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 149 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.149.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.149.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 150 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.150.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.150.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 151 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.151.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.151.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 152 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.152.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.152.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 153 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.153.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.153.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 154 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.154.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.154.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 155 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.155.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.155.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 156 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.156.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.156.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 157 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.157.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.157.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 158 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.158.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.158.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 159 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.159.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.159.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 160 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.160.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.160.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 161 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.161.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.161.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 162 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.162.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.162.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 163 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.163.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.163.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 164 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.164.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.164.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 165 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.165.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.165.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 166 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.166.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.166.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 167 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.167.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.167.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 168 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.168.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.168.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 169 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.169.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.169.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 170 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.170.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.170.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 171 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.171.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.171.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 172 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.172.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.172.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 173 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.173.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.173.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 174 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.174.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.174.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 175 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.175.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.175.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 176 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.176.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.176.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 177 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.177.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.177.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 178 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.178.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.178.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 179 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.179.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.179.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 180 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.180.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.180.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 181 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.181.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.181.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 182 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.182.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.182.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 183 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.183.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.183.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 184 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.184.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.184.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 185 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.185.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.185.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 186 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.186.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.186.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 187 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.187.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.187.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 188 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.188.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.188.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 189 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.189.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.189.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 190 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.190.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.190.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 191 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.191.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.191.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 192 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.192.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.192.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 193 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.193.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.193.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 194 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.194.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.194.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 195 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.195.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.195.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 196 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.196.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.196.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 197 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.197.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.197.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 198 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.198.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.198.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.198.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.198.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.198.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 199 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.199.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.199.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.199.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.199.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.199.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 200 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.200.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.200.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.200.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.200.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.200.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 201 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.201.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.201.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.201.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.201.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.201.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 202 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.202.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.202.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.202.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.202.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.202.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 203 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.203.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.203.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.203.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.203.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.203.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 204 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.204.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.204.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.204.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.204.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.204.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 205 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.205.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.205.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.205.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.205.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.205.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 206 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.206.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.206.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.206.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.206.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.206.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 207 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.207.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.207.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.207.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.207.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.207.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 208 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.208.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.208.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.208.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.208.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.208.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 209 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.209.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.209.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.209.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.209.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.209.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 210 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.210.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.210.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.210.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.210.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.210.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 211 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.211.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.211.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.211.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.211.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.211.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 212 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.212.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.212.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.212.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.212.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.212.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 213 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.213.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.213.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.213.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.213.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.213.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 214 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.214.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.214.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.214.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.214.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.214.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 215 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.215.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.215.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.215.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.215.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.215.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 216 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.216.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.216.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.216.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.216.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.216.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 217 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.217.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.217.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.217.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.217.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.217.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 218 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.218.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.218.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.218.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.218.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.218.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 219 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.219.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.219.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.219.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.219.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.219.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 220 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.220.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.220.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.220.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.220.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.220.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 221 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.221.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.221.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.221.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.221.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.221.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 222 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.222.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.222.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.222.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.222.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.222.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 223 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.223.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.223.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.223.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.223.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.223.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 224 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.224.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.224.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.224.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.224.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.224.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 225 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.225.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.225.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.225.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.225.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.225.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 226 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.226.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.226.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.226.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.226.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.226.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 227 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.227.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.227.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.227.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.227.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.227.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 228 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.228.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.228.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.228.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.228.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.228.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 229 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.229.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.229.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.229.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.229.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.229.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 230 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.230.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.230.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.230.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.230.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.230.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 231 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.231.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.231.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.231.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.231.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.231.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 232 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.232.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.232.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.232.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.232.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.232.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 233 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.233.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.233.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.233.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.233.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.233.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 234 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.234.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.234.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.234.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.234.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.234.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 235 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.235.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.235.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.235.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.235.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.235.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 236 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.236.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.236.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.236.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.236.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.236.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 237 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.237.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.237.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.237.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.237.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.237.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 238 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.238.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.238.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.238.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.238.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.238.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 239 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.239.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.239.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.239.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.239.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.239.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 240 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.240.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.240.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.240.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.240.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.240.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 241 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.241.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.241.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.241.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.241.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.241.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 242 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.242.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.242.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.242.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.242.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.242.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 243 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.243.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.243.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.243.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.243.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.243.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 244 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.244.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.244.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.244.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.244.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.244.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 245 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.245.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.245.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.245.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.245.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.245.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 246 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.246.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.246.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.246.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.246.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.246.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 247 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.247.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.247.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.247.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.247.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.247.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 248 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.248.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.248.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.248.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.248.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.248.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 249 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.249.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.249.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.249.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.249.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.249.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 250 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.250.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.250.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.250.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.250.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.250.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 251 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.251.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.251.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.251.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.251.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.251.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 252 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.252.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.252.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.252.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.252.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.252.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 253 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.253.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.253.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.253.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.253.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.253.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 254 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.254.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.254.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.254.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.254.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.254.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 255 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.255.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.255.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.255.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.255.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.255.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 256 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.256.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.256.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.256.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.256.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.256.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 257 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.257.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.257.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.257.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.257.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.257.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 258 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.258.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.258.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.258.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.258.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.258.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

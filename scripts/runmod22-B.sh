export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 100 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.100.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.100.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.100.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.100.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.100.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 101 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.101.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.101.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.101.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.101.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.101.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 102 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.102.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.102.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.102.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.102.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.102.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 103 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.103.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.103.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.103.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.103.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.103.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 104 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.104.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.104.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.104.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.104.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.104.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 105 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.105.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.105.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.105.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.105.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.105.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 106 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.106.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.106.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.106.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.106.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.106.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 107 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.107.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.107.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.107.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.107.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.107.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 108 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.108.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.108.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.108.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.108.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.108.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 109 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.109.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.109.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.109.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.109.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.109.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 110 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.110.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.110.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.110.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.110.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.110.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 111 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.111.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.111.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.111.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.111.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.111.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 112 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.112.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.112.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.112.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.112.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.112.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 113 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.113.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.113.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.113.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.113.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.113.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 114 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.114.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.114.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.114.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.114.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.114.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 115 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.115.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.115.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.115.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.115.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.115.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 116 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.116.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.116.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.116.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.116.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.116.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 117 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.117.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.117.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.117.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.117.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.117.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 118 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.118.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.118.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.118.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.118.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.118.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 119 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.119.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.119.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.119.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.119.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.119.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 120 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.120.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.120.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.120.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.120.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.120.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 121 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.121.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.121.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.121.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.121.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.121.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 122 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.122.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.122.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.122.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.122.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.122.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 123 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.123.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.123.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.123.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.123.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.123.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 124 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.124.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.124.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.124.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.124.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.124.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 125 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.125.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.125.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.125.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.125.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.125.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 126 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.126.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.126.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.126.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.126.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.126.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 127 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.127.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.127.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.127.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.127.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.127.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 128 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.128.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.128.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.128.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.128.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.128.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 129 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.129.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.129.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.129.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.129.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.129.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 130 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.130.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.130.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.130.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 131 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.131.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.131.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.131.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 132 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.132.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.132.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.132.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 133 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.133.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.133.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.133.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 134 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.134.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.134.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.134.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 135 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.135.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.135.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.135.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 136 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.136.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.136.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.136.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 137 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.137.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.137.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.137.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 138 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.138.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.138.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.138.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 139 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.139.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.139.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.139.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 140 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.140.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.140.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.140.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 141 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.141.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.141.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.141.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 142 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.142.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.142.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.142.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 143 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.143.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.143.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.143.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 144 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.144.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.144.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.144.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 145 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.145.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.145.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.145.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.145.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.145.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 146 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.146.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.146.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.146.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.146.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.146.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 147 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.147.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.147.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.147.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.147.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.147.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 148 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.148.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.148.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.148.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.148.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.148.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 149 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.149.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 150 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.150.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 151 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.151.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 152 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.152.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 153 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.153.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 154 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.154.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 155 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.155.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 156 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.156.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 157 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.157.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 158 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.158.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 159 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.159.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 160 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.160.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 161 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.161.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 162 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.162.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 163 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.163.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 164 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.164.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 165 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.165.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 166 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.166.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 167 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.167.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 168 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.168.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 169 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.169.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 170 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.170.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 171 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.171.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 172 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.172.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 173 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.173.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 174 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.174.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 175 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.175.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 176 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.176.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 177 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.177.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 178 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.178.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 179 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.179.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 180 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.180.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 181 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.181.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 182 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.182.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 183 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.183.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 184 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.184.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 185 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.185.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 186 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.186.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 187 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.187.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 188 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.188.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 189 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.189.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 190 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.190.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 191 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.191.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 192 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.192.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 193 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.193.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 194 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.194.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 195 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.195.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 196 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.196.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 197 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.197.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

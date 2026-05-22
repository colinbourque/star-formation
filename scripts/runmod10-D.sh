export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 388 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.388.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.388.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.388.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.388.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.388.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 389 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.389.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.389.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.389.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.389.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.389.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 390 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.390.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.390.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.390.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.390.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.390.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 391 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.391.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.391.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.391.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.391.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.391.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 392 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.392.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.392.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.392.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.392.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.392.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 393 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.393.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.393.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.393.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.393.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.393.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 394 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.394.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.394.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.394.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.394.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.394.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 395 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.395.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.395.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.395.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.395.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.395.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 396 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.396.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.396.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.396.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.396.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.396.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 397 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.397.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.397.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.397.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.397.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.397.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 398 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.398.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.398.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.398.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.398.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.398.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 399 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.399.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.399.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.399.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.399.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.399.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 400 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.400.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.400.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.400.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.400.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.400.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 401 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.401.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.401.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.401.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.401.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.401.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 402 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.402.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.402.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.402.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.402.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.402.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 403 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.403.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.403.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.403.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.403.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.403.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 404 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.404.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.404.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.404.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.404.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.404.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 405 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.405.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.405.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.405.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.405.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.405.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 406 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.406.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.406.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.406.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.406.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.406.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 407 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.407.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.407.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.407.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.407.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.407.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 408 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.408.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.408.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.408.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.408.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.408.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 409 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.409.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.409.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.409.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.409.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.409.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 410 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.410.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.410.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.410.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.410.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.410.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 411 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.411.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.411.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.411.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.411.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.411.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 412 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.412.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.412.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.412.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.412.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.412.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 413 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.413.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.413.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.413.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.413.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.413.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 414 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.414.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.414.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.414.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.414.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.414.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 415 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.415.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.415.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.415.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.415.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.415.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 416 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.416.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.416.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.416.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.416.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.416.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 417 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.417.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.417.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.417.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.417.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.417.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 418 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.418.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.418.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.418.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.418.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.418.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 419 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.419.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.419.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.419.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.419.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.419.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 420 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.420.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.420.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.420.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.420.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.420.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 421 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.421.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.421.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.421.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.421.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.421.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 422 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.422.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.422.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.422.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.422.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.422.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 423 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.423.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.423.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.423.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.423.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.423.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 424 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.424.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.424.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.424.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.424.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.424.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 425 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.425.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.425.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.425.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.425.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.425.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 426 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.426.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.426.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.426.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.426.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.426.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 427 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.427.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.427.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.427.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.427.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.427.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 428 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.428.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.428.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.428.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.428.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.428.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 429 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.429.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.429.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.429.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.429.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.429.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 430 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.430.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.430.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.430.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.430.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.430.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 431 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.431.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.431.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.431.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.431.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.431.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 432 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.432.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.432.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.432.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.432.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.432.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 433 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.433.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.433.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.433.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.433.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.433.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 434 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.434.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.434.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.434.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.434.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.434.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 435 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.435.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.435.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.435.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.435.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.435.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 436 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.436.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.436.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.436.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.436.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.436.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 437 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.437.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.437.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.437.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.437.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.437.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 438 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.438.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.438.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.438.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.438.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.438.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 439 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.439.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.439.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.439.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.439.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.439.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 440 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.440.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.440.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.440.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.440.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.440.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 441 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.441.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.441.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.441.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.441.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.441.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 442 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.442.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.442.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.442.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.442.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.442.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 443 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.443.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.443.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.443.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.443.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.443.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 444 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.444.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.444.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.444.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.444.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.444.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 445 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.445.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.445.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.445.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.445.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.445.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 446 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.446.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.446.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.446.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.446.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.446.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 447 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.447.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.447.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.447.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.447.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.447.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 448 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.448.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.448.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.448.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.448.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.448.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 449 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.449.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.449.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.449.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.449.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.449.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 450 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.450.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.450.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.450.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.450.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.450.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 451 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.451.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.451.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.451.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.451.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.451.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 452 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.452.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.452.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.452.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.452.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.452.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 453 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.453.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.453.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.453.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.453.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.453.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 454 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.454.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.454.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.454.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.454.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.454.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 455 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.455.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.455.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.455.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.455.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.455.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 456 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.456.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.456.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.456.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.456.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.456.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 457 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.457.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.457.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.457.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.457.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.457.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 458 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.458.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.458.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.458.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.458.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.458.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 459 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.459.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.459.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.459.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.459.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.459.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 460 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.460.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.460.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.460.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.460.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.460.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 461 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.461.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.461.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.461.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.461.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.461.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 462 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.462.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.462.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.462.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.462.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.462.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 463 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.463.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.463.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.463.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.463.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.463.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 464 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.464.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.464.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.464.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.464.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.464.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 465 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.465.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.465.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.465.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.465.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.465.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 466 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.466.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.466.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.466.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.466.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.466.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 467 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.467.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.467.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.467.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.467.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.467.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 468 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.468.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.468.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.468.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.468.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.468.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 469 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.469.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.469.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.469.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.469.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.469.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 470 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.470.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.470.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.470.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.470.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.470.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 471 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.471.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.471.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.471.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.471.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.471.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 472 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.472.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.472.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.472.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.472.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.472.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 473 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.473.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.473.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.473.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.473.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.473.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 474 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.474.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.474.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.474.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.474.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.474.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 475 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.475.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.475.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.475.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.475.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.475.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 476 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.476.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.476.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.476.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.476.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.476.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 477 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.477.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.477.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.477.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.477.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.477.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 478 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.478.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.478.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.478.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.478.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.478.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 479 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.479.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.479.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.479.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.479.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.479.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 480 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.480.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.480.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.480.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.480.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.480.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 481 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.481.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.481.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.481.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.481.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.481.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 482 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.482.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.482.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.482.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.482.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.482.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 483 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.483.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.483.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.483.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.483.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.483.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 484 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.484.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.484.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.484.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.484.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.484.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 485 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.485.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.485.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.485.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.485.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.485.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 486 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.486.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.486.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.486.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.486.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.486.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 487 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.487.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.487.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.487.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.487.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.487.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 488 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.488.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.488.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.488.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.488.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.488.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 489 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.489.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.489.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.489.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.489.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.489.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 490 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.490.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.490.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.490.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.490.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.490.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 491 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.491.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.491.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.491.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.491.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.491.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 492 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.492.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.492.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.492.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.492.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.492.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 493 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.493.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.493.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.493.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.493.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.493.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 494 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.494.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.494.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.494.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.494.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.494.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 495 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.495.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.495.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.495.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.495.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.495.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 496 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.496.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.496.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.496.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.496.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.496.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 497 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.497.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.497.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.497.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.497.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.497.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 498 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.498.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.498.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.498.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.498.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.498.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 499 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.499.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.499.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.499.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.499.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.499.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 500 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.500.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.500.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.500.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.500.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.500.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 501 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.501.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.501.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.501.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.501.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.501.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 502 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.502.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.502.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.502.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.502.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.502.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 503 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.503.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.503.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.503.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.503.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.503.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 504 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.504.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.504.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.504.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.504.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.504.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 505 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.505.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.505.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.505.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.505.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.505.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 506 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.506.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.506.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.506.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.506.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.506.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 507 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.507.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.507.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.507.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.507.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.507.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 508 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.508.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.508.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.508.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.508.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.508.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 509 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.509.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.509.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.509.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.509.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.509.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 510 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.510.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.510.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.510.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.510.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.510.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 511 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.511.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.511.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.511.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.511.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.511.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 512 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.512.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.512.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.512.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.512.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.512.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 513 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.513.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.513.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.513.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.513.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.513.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 514 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.514.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.514.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.514.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.514.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.514.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 515 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.515.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.515.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.515.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.515.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.515.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 516 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.516.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.516.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.516.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.516.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.516.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

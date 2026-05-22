export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 259 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.259.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.259.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.259.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.259.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.259.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 260 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.260.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.260.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.260.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.260.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.260.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 261 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.261.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.261.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.261.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.261.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.261.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 262 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.262.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.262.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.262.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.262.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.262.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 263 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.263.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.263.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.263.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.263.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.263.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 264 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.264.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.264.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.264.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.264.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.264.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 265 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.265.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.265.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.265.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.265.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.265.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 266 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.266.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.266.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.266.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.266.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.266.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 267 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.267.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.267.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.267.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.267.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.267.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 268 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.268.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.268.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.268.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.268.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.268.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 269 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.269.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.269.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.269.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.269.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.269.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 270 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.270.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.270.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.270.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.270.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.270.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 271 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.271.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.271.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.271.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.271.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.271.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 272 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.272.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.272.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.272.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.272.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.272.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 273 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.273.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.273.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.273.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.273.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.273.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 274 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.274.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.274.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.274.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.274.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.274.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 275 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.275.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.275.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.275.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.275.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.275.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 276 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.276.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.276.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.276.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.276.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.276.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 277 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.277.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.277.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.277.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.277.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.277.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 278 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.278.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.278.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.278.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.278.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.278.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 279 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.279.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.279.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.279.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.279.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.279.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 280 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.280.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.280.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.280.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.280.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.280.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 281 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.281.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.281.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.281.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.281.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.281.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 282 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.282.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.282.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.282.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.282.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.282.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 283 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.283.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.283.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.283.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.283.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.283.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 284 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.284.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.284.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.284.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.284.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.284.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 285 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.285.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.285.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.285.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.285.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.285.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 286 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.286.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.286.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.286.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.286.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.286.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 287 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.287.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.287.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.287.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.287.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.287.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 288 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.288.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.288.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.288.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.288.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.288.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 289 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.289.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.289.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.289.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.289.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.289.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 290 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.290.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.290.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.290.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.290.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.290.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 291 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.291.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.291.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.291.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.291.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.291.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 292 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.292.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.292.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.292.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.292.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.292.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 293 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.293.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.293.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.293.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.293.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.293.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 294 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.294.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.294.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.294.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.294.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.294.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 295 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.295.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.295.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.295.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.295.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.295.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 296 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.296.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.296.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.296.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.296.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.296.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 297 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.297.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.297.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.297.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.297.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.297.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 298 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.298.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.298.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.298.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.298.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.298.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 299 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.299.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.299.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.299.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.299.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.299.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 300 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.300.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.300.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.300.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.300.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.300.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 301 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.301.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.301.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.301.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.301.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.301.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 302 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.302.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.302.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.302.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.302.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.302.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 303 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.303.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.303.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.303.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.303.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.303.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 304 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.304.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.304.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.304.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.304.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.304.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 305 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.305.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.305.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.305.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.305.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.305.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 306 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.306.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.306.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.306.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.306.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.306.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 307 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.307.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.307.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.307.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.307.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.307.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 308 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.308.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.308.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.308.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.308.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.308.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 309 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.309.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.309.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.309.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.309.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.309.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 310 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.310.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.310.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.310.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.310.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.310.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 311 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.311.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.311.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.311.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.311.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.311.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 312 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.312.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.312.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.312.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.312.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.312.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 313 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.313.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.313.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.313.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.313.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.313.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 314 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.314.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.314.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.314.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.314.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.314.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 315 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.315.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.315.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.315.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.315.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.315.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 316 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.316.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.316.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.316.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.316.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.316.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 317 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.317.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.317.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.317.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.317.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.317.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 318 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.318.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.318.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.318.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.318.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.318.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 319 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.319.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.319.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.319.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.319.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.319.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 320 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.320.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.320.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.320.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.320.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.320.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 321 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.321.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.321.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.321.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.321.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.321.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 322 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.322.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.322.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.322.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.322.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.322.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 323 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.323.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.323.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.323.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.323.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.323.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 324 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.324.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.324.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.324.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.324.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.324.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 325 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.325.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.325.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.325.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.325.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.325.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 326 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.326.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.326.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.326.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.326.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.326.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 327 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.327.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.327.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.327.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.327.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.327.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 328 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.328.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.328.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.328.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.328.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.328.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 329 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.329.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.329.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.329.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.329.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.329.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 330 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.330.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.330.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.330.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.330.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.330.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 331 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.331.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.331.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.331.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.331.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.331.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 332 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.332.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.332.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.332.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.332.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.332.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 333 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.333.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.333.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.333.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.333.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.333.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 334 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.334.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.334.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.334.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.334.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.334.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 335 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.335.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.335.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.335.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.335.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.335.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 336 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.336.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.336.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.336.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.336.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.336.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 337 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.337.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.337.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.337.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.337.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.337.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 338 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.338.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.338.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.338.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.338.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.338.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 339 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.339.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.339.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.339.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.339.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.339.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 340 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.340.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.340.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.340.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.340.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.340.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 341 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.341.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.341.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.341.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.341.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.341.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 342 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.342.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.342.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.342.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.342.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.342.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 343 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.343.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.343.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.343.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.343.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.343.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 344 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.344.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.344.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.344.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.344.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.344.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 345 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.345.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.345.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.345.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.345.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.345.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 346 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.346.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.346.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.346.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.346.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.346.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 347 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.347.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.347.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.347.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.347.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.347.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 348 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.348.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.348.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.348.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.348.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.348.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 349 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.349.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.349.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.349.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.349.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.349.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 350 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.350.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.350.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.350.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.350.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.350.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 351 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.351.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.351.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.351.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.351.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.351.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 352 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.352.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.352.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.352.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.352.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.352.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 353 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.353.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.353.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.353.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.353.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.353.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 354 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.354.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.354.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.354.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.354.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.354.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 355 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.355.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.355.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.355.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.355.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.355.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 356 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.356.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.356.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.356.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.356.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.356.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 357 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.357.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.357.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.357.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.357.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.357.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 358 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.358.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.358.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.358.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.358.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.358.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 359 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.359.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.359.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.359.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.359.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.359.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 360 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.360.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.360.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.360.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.360.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.360.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 361 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.361.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.361.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.361.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.361.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.361.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 362 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.362.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.362.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.362.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.362.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.362.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 363 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.363.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.363.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.363.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.363.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.363.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 364 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.364.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.364.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.364.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.364.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.364.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 365 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.365.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.365.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.365.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.365.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.365.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 366 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.366.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.366.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.366.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.366.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.366.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 367 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.367.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.367.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.367.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.367.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.367.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 368 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.368.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.368.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.368.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.368.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.368.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 369 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.369.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.369.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.369.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.369.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.369.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 370 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.370.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.370.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.370.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.370.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.370.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 371 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.371.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.371.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.371.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.371.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.371.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 372 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.372.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.372.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.372.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.372.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.372.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 373 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.373.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.373.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.373.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.373.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.373.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 374 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.374.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.374.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.374.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.374.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.374.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 375 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.375.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.375.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.375.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.375.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.375.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 376 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.376.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.376.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.376.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.376.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.376.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 377 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.377.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.377.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.377.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.377.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.377.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 378 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.378.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.378.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.378.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.378.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.378.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 379 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.379.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.379.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.379.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.379.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.379.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 380 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.380.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.380.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.380.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.380.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.380.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 381 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.381.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.381.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.381.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.381.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.381.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 382 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.382.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.382.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.382.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.382.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.382.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 383 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.383.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.383.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.383.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.383.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.383.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 384 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.384.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.384.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.384.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.384.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.384.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 385 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.385.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.385.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.385.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.385.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.385.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 386 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.386.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.386.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.386.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.386.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.386.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 387 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/stars.387.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/amr_grid.387.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/INPUT/dust_density.387.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/dust_temperature.387.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model10/OUTPUT/spectrum.387.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

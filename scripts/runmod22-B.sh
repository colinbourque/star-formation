export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 75 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.075.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.075.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.075.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.075.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.075.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 76 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.076.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.076.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.076.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.076.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.076.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 77 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.077.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.077.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.077.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.077.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.077.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 78 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.078.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.078.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.078.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.078.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.078.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 79 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.079.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.079.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.079.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.079.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.079.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 80 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.080.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.080.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.080.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.080.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.080.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 81 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.081.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.081.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.081.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.081.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.081.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 82 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.082.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.082.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.082.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.082.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.082.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 83 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.083.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.083.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.083.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.083.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.083.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 84 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.084.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.084.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.084.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.084.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.084.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 85 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.085.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.085.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.085.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.085.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.085.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 86 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.086.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.086.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.086.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.086.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.086.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 87 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.087.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.087.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.087.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.087.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.087.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 88 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.088.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.088.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.088.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.088.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.088.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 89 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.089.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.089.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.089.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.089.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.089.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 90 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.090.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.090.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.090.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.090.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.090.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 91 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.091.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.091.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.091.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.091.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.091.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 92 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.092.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.092.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.092.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.092.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.092.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 93 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.093.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.093.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.093.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.093.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.093.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 94 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.094.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.094.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.094.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.094.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.094.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 95 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.095.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.095.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.095.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.095.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.095.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 96 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.096.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.096.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.096.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.096.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.096.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 97 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.097.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.097.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.097.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.097.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.097.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 98 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.098.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.098.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.098.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.098.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.098.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 99 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.099.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.099.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.099.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.099.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-B/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.099.incl85.out 
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
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

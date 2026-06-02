export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 149 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.149.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.149.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.149.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 150 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.150.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.150.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.150.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 151 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.151.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.151.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.151.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 152 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.152.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.152.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.152.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 153 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.153.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.153.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.153.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 154 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.154.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.154.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.154.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 155 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.155.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.155.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.155.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 156 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.156.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.156.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.156.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 157 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.157.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.157.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.157.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 158 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.158.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.158.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.158.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 159 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.159.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.159.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.159.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 160 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.160.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.160.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.160.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 161 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.161.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.161.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.161.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 162 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.162.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.162.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.162.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 163 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.163.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.163.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.163.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 164 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.164.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.164.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.164.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 165 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.165.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.165.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.165.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 166 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.166.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.166.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.166.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 167 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.167.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.167.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.167.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 168 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.168.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.168.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.168.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 169 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.169.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.169.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.169.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 170 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.170.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.170.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.170.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 171 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.171.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.171.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.171.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 172 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.172.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.172.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.172.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 173 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.173.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.173.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.173.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 174 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.174.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.174.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.174.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 175 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.175.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.175.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.175.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 176 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.176.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.176.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.176.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 177 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.177.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.177.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.177.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 178 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.178.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.178.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.178.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 179 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.179.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.179.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.179.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 180 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.180.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.180.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.180.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 181 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.181.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.181.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.181.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 182 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.182.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.182.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.182.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 183 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.183.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.183.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.183.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 184 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.184.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.184.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.184.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 185 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.185.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.185.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.185.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 186 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.186.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.186.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.186.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 187 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.187.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.187.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.187.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 188 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.188.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.188.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.188.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 189 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.189.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.189.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.189.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 190 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.190.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.190.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.190.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 191 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.191.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.191.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.191.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 192 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.192.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.192.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.192.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 193 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.193.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.193.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.193.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 194 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.194.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.194.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.194.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 195 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.195.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.195.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.195.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 196 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.196.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.196.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.196.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 197 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.197.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.197.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.197.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 198 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.198.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.198.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.198.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.198.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.198.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 199 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.199.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.199.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.199.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.199.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.199.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 200 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.200.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.200.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.200.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.200.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.200.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 201 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.201.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.201.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.201.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.201.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.201.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 202 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.202.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.202.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.202.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.202.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.202.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 203 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.203.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.203.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.203.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.203.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.203.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 204 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.204.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.204.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.204.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.204.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.204.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 205 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.205.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.205.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.205.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.205.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.205.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 206 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.206.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.206.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.206.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.206.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.206.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 207 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.207.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.207.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.207.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.207.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.207.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 208 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.208.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.208.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.208.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.208.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.208.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 209 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.209.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.209.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.209.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.209.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.209.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 210 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.210.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.210.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.210.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.210.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.210.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 211 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.211.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.211.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.211.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.211.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.211.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 212 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.212.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.212.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.212.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.212.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.212.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 213 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.213.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.213.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.213.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.213.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.213.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 214 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.214.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.214.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.214.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.214.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.214.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 215 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.215.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.215.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.215.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.215.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.215.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 216 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.216.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.216.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.216.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.216.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.216.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 217 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.217.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.217.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.217.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.217.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.217.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 218 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.218.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.218.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.218.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.218.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.218.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 219 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.219.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.219.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.219.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.219.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.219.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 220 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.220.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.220.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.220.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.220.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.220.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 221 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.221.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.221.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.221.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.221.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.221.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 222 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.222.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.222.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.222.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.222.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-C/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.222.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

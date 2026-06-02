export STARTTIME=$(date) 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/aperture_info.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/aperture_info.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustopac.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dustopac.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/radmc3d.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/radmc3d.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/wavelength_micron.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/wavelength_micron.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/dustkappa_OH5.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dustkappa_OH5.inp 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 223 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.223.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.223.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.223.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.223.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.223.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 224 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.224.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.224.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.224.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.224.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.224.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 225 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.225.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.225.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.225.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.225.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.225.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 226 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.226.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.226.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.226.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.226.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.226.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 227 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.227.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.227.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.227.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.227.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.227.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 228 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.228.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.228.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.228.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.228.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.228.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 229 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.229.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.229.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.229.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.229.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.229.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 230 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.230.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.230.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.230.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.230.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.230.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 231 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.231.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.231.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.231.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.231.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.231.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 232 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.232.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.232.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.232.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.232.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.232.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 233 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.233.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.233.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.233.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.233.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.233.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 234 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.234.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.234.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.234.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.234.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.234.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 235 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.235.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.235.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.235.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.235.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.235.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 236 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.236.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.236.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.236.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.236.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.236.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 237 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.237.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.237.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.237.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.237.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.237.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 238 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.238.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.238.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.238.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.238.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.238.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 239 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.239.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.239.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.239.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.239.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.239.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 240 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.240.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.240.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.240.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.240.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.240.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 241 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.241.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.241.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.241.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.241.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.241.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 242 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.242.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.242.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.242.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.242.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.242.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 243 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.243.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.243.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.243.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.243.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.243.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 244 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.244.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.244.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.244.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.244.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.244.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 245 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.245.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.245.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.245.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.245.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.245.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 246 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.246.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.246.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.246.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.246.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.246.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 247 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.247.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.247.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.247.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.247.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.247.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 248 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.248.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.248.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.248.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.248.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.248.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 249 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.249.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.249.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.249.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.249.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.249.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 250 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.250.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.250.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.250.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.250.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.250.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 251 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.251.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.251.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.251.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.251.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.251.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 252 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.252.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.252.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.252.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.252.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.252.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 253 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.253.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.253.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.253.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.253.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.253.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 254 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.254.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.254.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.254.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.254.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.254.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 255 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.255.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.255.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.255.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.255.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.255.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 256 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.256.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.256.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.256.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.256.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.256.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 257 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.257.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.257.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.257.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.257.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.257.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 258 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.258.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.258.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.258.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.258.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.258.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 259 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.259.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.259.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.259.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.259.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.259.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 260 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.260.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.260.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.260.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.260.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.260.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 261 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.261.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.261.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.261.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.261.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.261.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 262 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.262.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.262.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.262.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.262.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.262.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 263 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.263.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.263.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.263.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.263.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.263.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 264 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.264.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.264.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.264.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.264.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.264.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 265 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.265.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.265.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.265.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.265.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.265.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 266 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.266.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.266.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.266.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.266.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.266.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 267 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.267.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.267.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.267.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.267.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.267.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 268 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.268.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.268.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.268.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.268.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.268.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 269 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.269.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.269.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.269.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.269.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.269.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 270 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.270.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.270.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.270.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.270.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.270.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 271 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.271.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.271.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.271.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.271.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.271.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 272 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.272.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.272.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.272.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.272.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.272.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 273 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.273.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.273.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.273.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.273.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.273.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 274 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.274.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.274.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.274.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.274.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.274.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 275 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.275.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.275.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.275.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.275.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.275.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 276 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.276.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.276.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.276.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.276.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.276.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 277 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.277.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.277.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.277.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.277.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.277.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 278 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.278.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.278.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.278.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.278.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.278.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 279 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.279.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.279.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.279.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.279.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.279.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 280 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.280.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.280.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.280.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.280.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.280.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 281 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.281.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.281.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.281.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.281.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.281.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 282 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.282.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.282.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.282.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.282.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.282.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 283 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.283.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.283.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.283.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.283.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.283.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 284 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.284.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.284.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.284.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.284.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.284.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 285 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.285.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.285.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.285.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.285.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.285.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 286 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.286.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.286.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.286.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.286.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.286.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 287 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.287.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.287.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.287.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.287.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.287.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 288 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.288.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.288.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.288.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.288.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.288.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 289 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.289.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.289.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.289.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.289.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.289.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 290 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.290.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.290.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.290.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.290.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.290.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 291 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.291.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.291.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.291.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.291.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.291.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 292 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.292.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.292.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.292.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.292.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.292.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 293 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.293.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.293.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.293.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.293.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.293.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 294 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.294.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.294.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.294.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.294.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.294.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 295 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.295.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.295.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.295.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.295.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.295.incl85.out 
echo --------------------------------------------------------- 
echo BEGIN MODEL STEP 296 
echo --------------------------------------------------------- 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/stars.296.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/stars.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/amr_grid.296.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/amr_grid.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/INPUT/dust_density.296.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_density.inp 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/AUX/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp 
radmc3d mctherm setthreads 4 countwrite 1000000 
cp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/dust_temperature.dat /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/dust_temperature.296.dat 
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/external_source.inp /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/no_external_source.inp 
radmc3d sed incl 5 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl05.out 
radmc3d sed incl 15 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl15.out 
radmc3d sed incl 25 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl25.out 
radmc3d sed incl 35 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl35.out 
radmc3d sed incl 45 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl45.out 
radmc3d sed incl 55 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl55.out 
radmc3d sed incl 65 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl65.out 
radmc3d sed incl 75 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl75.out 
radmc3d sed incl 85 useapert dpc 140 setthreads 4
mv /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/RUNDIR-D/spectrum.out /Users/cbourque/astrophysics/star-formation/SIMULATIONS/model22/OUTPUT/spectrum.296.incl85.out 
echo STARTED AT $STARTTIME 
echo ENDED AT $(date) 

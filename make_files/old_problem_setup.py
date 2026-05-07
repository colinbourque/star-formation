import numpy as np
from radmc_utils import *

au      = 1.49598e13     # Astronomical Unit       [cm]
gg      = 6.67430e-8     ## gravitational constant, G, in cgs units (dyn cm^2 g^-2)
Delta_Q = -3.52e-4       # Relates x and y based on how far from spherical symmetry
ms      = 1.98892e33     # Solar mass              [g]
ts      = 5.78e3         # Solar temperature       [K]
ls      = 3.8525e33      # Solar luminosity        [erg/s]
rs      = 6.96e10        # Solar radius            [cm]
pc      = 3.08567758e18  # cm per pc according to astropy.units
c       = 299792458*100  # speed of light in cm/s
ssb     = 5.670374419e-5 # stefan-boltzmann constant [erg⋅cm−2⋅s−1⋅K−4] (this does not get used after all)

simdir = '/Users/cbourque/astrophysics/star-formation/SIMULATIONS'

get_modelnum = True
while get_modelnum:
  modelnumber = input('Generate input files for model number...')
  try:
    modelnumber = int(modelnumber)
    get_modelnum = False
  except:
    print('Input not an integer')
print(f'Generating input files for model{modelnumber:02n}')

## write the dust opacity control file
with open(f'{simdir}/model{modelnumber:02n}/AUX/dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('OH5             Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')

## write the radmc3d.inp control file
nphot = 10000000
with open(f'{simdir}/model{modelnumber:02n}/AUX/radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 1\n')   # Put this to 1 for isotropic scattering
    f.write('istar_sphere = 0\n')          ## point-source stars matches 2d code
    f.write('iranfreqmode = 1\n')
    f.write('modified_random_walk = 1\n')  ## this analytically determines outward travel from high-τ cells

## Get wavelength information (for writing the stars.inp file)
nlam=100
with open(f'aux/wavelength_micron.inp','r') as f:
    wavs = f.readlines()
wavs = [float(wav) for wav in wavs[1:]]

## Produce a control file defining the wavelength grid
with open(f'{simdir}/model{modelnumber:02n}/AUX/wavelength_micron.inp', 'w+') as w:
    w.write('100\n')
    for wav in wavs:
        w.write(f'{wav:.6e}\n')

## Produce a control file specifying the absorption and scattering opacities of OH5 dust
with open('dustopac_1.inp') as f:
    opacs = f.readlines()
absp = opacs[2:102]
scat = opacs[102:]
scat = [float(cap) for cap in scat]
absp = [float(cap) for cap in absp]
with open(f'{simdir}/model{modelnumber:02n}/AUX/dustkappa_OH5.inp', 'w+') as f:
    f.write(f'2\n')
    f.write(f'100\n')
    for i in range(len(wavs)):
        f.write(f'{wavs[i]:.6e}\t{absp[i]:.6e}\t{scat[i]:.6e}\n')

## Produce a control file specifying an external radiation field without the CMB long-wavelength peak
ext_int = np.loadtxt('external_meanint.inp', skiprows=2).T
intens = ext_int[1]
zero = 0
with open(f'{simdir}/model{modelnumber:02n}/AUX/external_source.inp', 'w+') as f:
    f.write(f'2\n')
    f.write(f'100\n')
    for i in range(len(wavs)):
        f.write(f'{wavs[i]:.6e}\n')
    for i in range(len(wavs)):
        f.write(f'{intens[i]:.6e}\n')
        # if wavs[i] >= 400:
        #     f.write(f'{zero:.6e}\n')
        # elif wavs[i] < 400:
        #     f.write(f'{intens[i]:.6e}\n')


## star is positioned at origin
pstar = [0,0,0]

aps = [2, 6, 20, 50]
## The following code has been commented out to match the actual apertures applied in Dunham & Voroboyov 2012
## It can be uncommented in order to apply a more physically-motivated long-wavelength aperture
# outermost_radius = 4125.2671*au
# dobserver = 140*pc
# vangle_rad = (outermost_radius*1.10)/dobserver
# vangle_arcsec = round(vangle_rad*206265, 1)
# aps = [2, 6, 20, vangle_arcsec]

with open(f'{simdir}/model{modelnumber:02n}/AUX/aperture_info.inp', 'w+') as f:
    f.write(f'1\n')
    f.write(f'100\n')
    for i in range(len(wavs)):
        if wavs[i] <= 10:
            f.write(f'{wavs[i]:.6e}\t{aps[0]:.6e}\n')
        elif 10 < wavs[i] <= 40:
            f.write(f'{wavs[i]:.6e}\t{aps[1]:.6e}\n')
        elif 40 < wavs[i] <= 100:
            f.write(f'{wavs[i]:.6e}\t{aps[2]:.6e}\n')
        elif 100 < wavs[i]:
            f.write(f'{wavs[i]:.6e}\t{aps[3]:.6e}\n')

## files input to the radmc3d runs will be written here
indir = f'{simdir}/model{modelnumber:02n}/INPUT'

#### Information which would go into each individual run
nr = 130
nt = 100
ntex = 20

## 2-dimensional simulations will always take a singular 0 < phi < 2π argument
nphi = 1
phii = np.array([0.0, np.pi*2])
# phir = np.pi/2.0 - tt

## Accretion disk parameters which may or may not ever change?
H0 = 10*au
s0 = 100*au
beta = 1.25
alpha = 2.25
rho0 = 1e-16

#### CREATE INPUT DATA
mpars = np.loadtxt(f'{simdir}/model_parameters/model_parameters_{modelnumber}.tbl', skiprows=1) 
mpars = mpars[~np.isnan(mpars).any(axis=1)].T ## remove entries with nan values

v_collapse = 0 
r_out = np.max([mpars[7][0], mpars[10][0]])*au ## starts the outer radius as being the greater of the disk or core

for i in range(1,len(mpars[0])):
    ## get specific values for this model step
    time_myr  = mpars[0][i]
    time_mt0  = mpars[1][i]
    mstar     = mpars[2][i]
    lstar     = mpars[3][i]
    rstar     = mpars[4][i]
    tstar     = mpars[5][i]
    rdisk_in  = mpars[6][i]
    rdisk_out = mpars[7][i]
    mdisk     = mpars[8][i]
    renv_in   = mpars[9][i]
    renv_out  = mpars[10][i]
    menv      = mpars[11][i]
    Omega_0   = mpars[12][i]
    c_s       = mpars[13][i]
    
    ## Set the outer boundary
    delta_t_myr = mpars[0][i] - mpars[0][i-1]
    if v_collapse > 0:
        delta_t = (delta_t_myr)*(1000000*365.*24.*60.*60.) # seconds
        r_out -= v_collapse*delta_t
    
    ## Set the inner boundary
    if rdisk_out == 0:   ## if no disk then only core
        r_in = renv_in*au
    if rdisk_out != 0:   ## if yes disk then disk or core inner, whichever is closer
        r_in = np.min([rdisk_in, renv_in])*au
    if rstar*rs >= r_in: ## if r_in within star... no it isn't
        r_in = rstar*rs*1.10
    
    if rdisk_out > r_out: 
        print('Disk outer radius is greater than the space being simulated')
        print('Colin does not know why this check is needed so consider him baffled if you see this')
        break
    
    ## Make the simulation mesh grids
    r_walls, r_means         = make_r_grid(nr=nr, rin=r_in, rout=r_out)
    theta_walls, theta_means = make_theta_grid(nt=nt, ntex=ntex, tmps=.5)
    rr, tt = np.meshgrid(r_means, theta_means, indexing='ij')
    phir = np.pi/2.0 - tt
    
    ## Get the volumes of each 'cell' for use in normalizing densities
    rl  = r_integ(r_walls[:-1], r_walls[1:])            ## Integrate "length" between cell walls
    thl = th_integ(theta_walls[:-1], theta_walls[1:])   ## "
    rrl, ththl = np.meshgrid(rl, thl, indexing='ij')    ## Mesh the "length" grids together
    cv = 2*np.pi*rrl*ththl                              ## 'cell volume' of every cell in the simulation-space

    ## Sanity check that we have preserved the desired dimensionality
    if not rr.shape == (nr,nt):
        print(f'Simulation grid does not match requested sizes:')
        print(f'Expected r x theta = {nr, nt} got {rr.shape}')
        break
        
    ## Make the accretion disk
    if not mdisk == 0:
        zz = rr*np.cos(tt)
        ss = np.sqrt(rr**2 - zz**2)
        Hs = H0*(ss/s0)**beta
        rhod = rho0 * (ss/s0)**(-alpha) * np.exp(-0.5*(zz/Hs)**2)

        rhod   = np.where(ss <= rdisk_out*au, rhod, 0) ## if not within disk bounds, rho -> 0

        temp_disk_mass = np.sum(rhod*cv)/ms   ## get current mass of disk (meaningless)
        rhod *= mdisk / temp_disk_mass  ## rescale to match hydro-determined disk mass
    
    elif mdisk == 0: 
        rhod = np.full((nr, nt), 0)
    
    
    ## Now make the core density profile
    rhoenv, u_r, y_tsc, x_tsc, tau_tsc = make_cdp(rr=rr, tt=tt, Omega_0=Omega_0, 
                                                  modeltime=time_mt0, Delta_Q=Delta_Q, c_s=c_s)
    
    rhoenv = np.where(rr >= renv_in*au, rhoenv, 0) ## if not within core bounds, rho -> 0
    
    temp_core_mass = np.sum(rhoenv*cv)/ms ## get current mass of core (meaningless)
    rhoenv *= menv / temp_core_mass       ## rescale to match hydro-determined core mass
    
    ##### reset v_collapse based on u_r IFF the collapse-front has reached the outer bound
    if u_r[nr-1,int(nt/2)] > 0:
        v_collapse = u_r[nr-1,int(nt/2)]
    
    ## Combine the density profiles in their respective regions
    # rhod   = np.where(rr <= rdisk_out*au, rhod, 0) ## if not within disk bounds, rho -> 0
    # rhoenv = np.where(rr >= renv_in*au, rhoenv, 0) ## if not within core bounds, rho -> 0
    rho = np.maximum(rhoenv, rhod)                      ## cell takes whichever density is greater
    
    rho *= (1/100)
    
    #### Determine parameters necessary for writing the radmc3d input files
    ## things like the star position would go here?
    ## although I think most of this ends up being hardcoded anyways
    
    #### Begin writing the necessary radmc3d input files
        
    ## write the stars.inp file
    pstar = [0,0,0]
    with open(f'{indir}/stars.{i:03n}.inp','w+') as f:
        f.write('2\n')
        f.write('1 %d\n\n'%(nlam))
        f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar*rs,mstar*ms,pstar[0],pstar[1],pstar[2]))
        for value in wavs:
            f.write('%13.6e\n'%(value))
        f.write('\n%13.6e\n'%(-tstar))

    ## write the grid file
    with open(f'{indir}/amr_grid.{i:03n}.inp','w+') as f:
        f.write('1\n')                       # iformat
        f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
        f.write('100\n')                     # Coordinate system
        f.write('0\n')                       # gridinfo
        f.write('1 1 0\n')                   # Include x,y,z coordinate
        f.write('%d %d %d\n'%(nr,nt,nphi))   # Size of grid
        for value in r_walls:
            f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
        for value in theta_walls:
            f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
        for value in phii:
            f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)

    ## write the density file
    with open(f'{indir}/dust_density.{i:03n}.inp','w+') as f:
        f.write('1\n')                       # Format number
        f.write('%d\n'%(nr*nt*nphi))           # Nr of cells
        f.write('1\n')                       # Nr of dust species
        data = rho.ravel(order='F')         # Create a 1-D view, fortran-style indexing
        data.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')

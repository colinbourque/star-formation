import numpy as np
au  = 1.49598e13     # Astronomical Unit       [cm]
gg  = 6.67430e-8     ## gravitational constant, G, in cgs units (dyn cm^2 g^-2)
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
pc  = 3.08567758e18  # cm per pc according to astropy.units
ssb = 5.670374419e-5 # stefan-boltzmann constant [erg⋅cm−2⋅s−1⋅K−4]
cc  = 299792458*100  # speed of light in cm/s



def make_r_grid(nr=130, rin=10*au, rout=100*au):
    ri = np.geomspace(rin, rout, nr+1)    ## list of cell boundaries
    rc = np.sqrt( ri[:-1]*ri[1:] )        ## geometric means for within each cell (rm == Radial Means)
    rc[0] = ri[0]                         ## redefine the innermost point to be equal to the inner boundary
    return ri, rc

def make_theta_grid(nt=40, ntex=10, tmps=.5):
    if ntex%2 == 1:
        print("Number of extra points not divisible by two")
        print("This may result in weird behavior")

    if (nt%1 != 0) or (ntex%1 != 0):
        print("Non-integer number of points requested. Rounding will occur.")

    nt   = int(nt)
    ntex = int(ntex) ## Need to sanitize these for use with integer division operator

    ti = np.array([*np.linspace(0,            np.pi/2-tmps, ntex//2, endpoint=False),
                   *np.linspace(np.pi/2-tmps, np.pi/2+tmps, nt-ntex+1),
                   *np.linspace(np.pi,        np.pi/2+tmps, ntex//2, endpoint=False)[::-1]])
    tc = 0.5*(ti[:-1]+ti[1:])
    return ti, tc

def cubesolve(b):
    x = np.full(4, -999.0)
    #pi = 3.14159265358979323846264338327950288
    pi = np.pi

    if abs(b[0]) < 1.0e-10:
        print("Coefficient of cubic term in CUBESOLVE < 1e-10")
        print("Aborting!!!")
        return np.full(4, np.nan)

    ## some parameters of probable numerical significance
    p=b[1]/b[0]
    q=b[2]/b[0]
    r=b[3]/b[0]
    ba=q-(p*p/3.)
    bb=((2.*(p**3))-(9.*p*q)+(27.*r))/27.0
    det=((bb*bb)/4.)+((ba**3)/27.)

    # find solution, check sign of 'det' to find root types

    # one real, two complex roots
    if det > 0.0:
        ca=((-bb/2.)+np.sqrt(det))
        if ca >= 0.0:
            sign=1.0
        if ca < 0.0:
            sign=-1.0
        capa=1.0*sign*(abs(ca)**(1./3.))
        cb=((-bb/2.)-np.sqrt(det))
        if cb >= 0.0:
            sign=1.0
        if cb < 0.0:
            sign=-1.0
        capb=1.0*sign*(abs(cb)**(1./3.))
        root1=capa+capb
        x[0]=root1-(p/3.)
        x[3]=1.0

    # three real distinct roots
    if det < 0.0:
        rg = (-bb/2.)/np.sqrt(-ba**(3./27.))
        phi = np.arccos(rg)
        cons = 2.0*np.sqrt(-ba/3.)
        x[0]=cons*np.cos(phi/3.)-(p/3.)
        x[1]=cons*np.cos((phi/3.)+(2.*pi/3.))-(p/3.)
        x[2]=cons*np.cos((phi/3.)+(4.*pi/3.))-(p/3.)
        x[3]=-1.0

    # three real roots
    if det == 0.0:
        bbtemp=-1.0*bb
        if bbtemp > 0.0:
            sign=1.0
        if bbtemp < 0.0:
            sign=-1.0
        capa=1.0*sign*(abs(-bb/2.)**(1./3.))
        x[0]=2.*capa-(p/3.)
        x[1]=(-1.0*capa)-(p/3.)
        x[3]=x[2]
        x[4]=0.0

    return x

def extrap(y_sim, y_tab, p_tab):
    y1, y2 = np.log10(np.array([y_tab[-2], y_tab[-1]]))
    p1, p2 = np.log10(np.array([p_tab[-2], p_tab[-1]]))
    m      = (p2-p1)/(y2-y1)
    pex    = m*(np.log10(y_sim)-y2)+p2
    return 10**pex

## Define two functions representing the indefinite integrals of each respective volume component,
## accepting the bounds as arguments
def r_integ(rl, rh):
    return -(1/3) * (rh**3 - rl**3)

def th_integ(thl, thh):
    return np.cos(thh) - np.cos(thl)


def make_cdp(rr, tt, Omega_0, modeltime, Delta_Q=-3.52e-4, c_s = 18800.0):

    ######################### read_dimensionless_grid.ipynb
    dd = np.loadtxt('grid.plt34.mod', skiprows=1).T ## transpose to fix indexing
    y_dd      = dd[0]
    alpha0_dd = dd[1]
    alphaM_dd = dd[2]
    alphaQ_dd = (4/3)*dd[3]
    V0_dd     = dd[4]
    VM_dd     = dd[5]
    VQ_dd     = (2/3)*dd[6]

    nr, nt = np.shape(rr)

    x_tsc   = rr/(c_s*modeltime*365.0*24.0*3600.0)
    tau_tsc = np.full((nr,nt), Omega_0*(modeltime*365.0*24.0*3600.0))
    p2      = (1./2.)*((3.0*(np.cos(tt)**2))-1.0)
    y_tsc   = x_tsc * (1.+((tau_tsc**2.0)*Delta_Q*p2))

    alpha_0_tsc = np.zeros((nr,nt))
    alpha_M_tsc = np.zeros((nr,nt))
    alpha_Q_tsc = np.zeros((nr,nt))
    V_0_tsc     = np.zeros((nr,nt))
    V_M_tsc     = np.zeros((nr,nt))
    V_Q_tsc     = np.zeros((nr,nt))
    rhoenv      = np.zeros((nr,nt))

    for i in range(nr):
        for j in range(nt):
            if y_tsc[i,j] < y_dd[6]:
                alpha_0_tsc[i,j] =  0
                alpha_M_tsc[i,j] =  0
                alpha_Q_tsc[i,j] =  0
                V_0_tsc[i,j] =  0
                V_M_tsc[i,j] =  0
                V_Q_tsc[i,j] =  0

            elif y_dd[6] <= y_tsc[i,j] <= y_dd[743]:
                alpha_0_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(alpha0_dd[6:744]))
                alpha_M_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(alphaM_dd[6:744]))
                alpha_Q_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(alphaQ_dd[6:744]))
                V_0_tsc[i,j] =  -10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(-V0_dd[6:744]))
                V_M_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(VM_dd[6:744]))
                V_Q_tsc[i,j] =  -10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(VQ_dd[6:744]))

            elif y_dd[743] < y_tsc[i,j] < y_dd[744]:
                alpha_0_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(alpha0_dd[6:744]))
                alpha_M_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(alphaM_dd[6:744]))
                alpha_Q_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(alphaQ_dd[6:744]))
                V_0_tsc[i,j] =  -10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(-V0_dd[6:744]))
                V_M_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[6:744]), np.log10(VM_dd[6:744]))
                V_Q_tsc[i,j] =  0

            elif y_dd[743] <= y_tsc[i,j] <= 1.:
                alpha_0_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(alpha0_dd[744:]))
                alpha_M_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(alphaM_dd[744:]))
                alpha_Q_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(alphaQ_dd[744:]))
                V_0_tsc[i,j] =  -10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:764]), np.log10(-V0_dd[744:764]))
                V_M_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:764]), np.log10(VM_dd[744:764]))
                V_Q_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(-VQ_dd[744:]))

            elif 1. < y_tsc[i,j] <= np.max(y_dd):
                alpha_0_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(alpha0_dd[744:]))
                alpha_M_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(alphaM_dd[744:]))
                alpha_Q_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(alphaQ_dd[744:]))
                V_0_tsc[i,j] =  0
                V_M_tsc[i,j] =  0
                V_Q_tsc[i,j] =  10**np.interp(np.log10(y_tsc[i,j]), np.log10(y_dd[744:]), np.log10(-VQ_dd[744:]))

            elif y_tsc[i,j] > np.max(y_dd):
                alpha_0_tsc[i,j] = extrap(y_tsc[i,j], y_dd[744:], alpha0_dd[744:])
                alpha_M_tsc[i,j] = extrap(y_tsc[i,j], y_dd[744:], alphaM_dd[744:])
                alpha_Q_tsc[i,j] = extrap(y_tsc[i,j], y_dd[744:], alphaQ_dd[744:])
                V_0_tsc[i,j] =  0
                V_M_tsc[i,j] =  0
                V_Q_tsc[i,j] = 1.*extrap(y_tsc[i,j], y_dd[744:], -VQ_dd[744:])

    u_r = c_s*( V_0_tsc + ( (tau_tsc**2.) * (V_M_tsc + (V_Q_tsc*p2) ) ) )

    for i in range(nr):
        for j in range(nt):

            # inner solution
            if x_tsc[i,j] < tau_tsc[i,j]**2:
                xmo=0.973546863
                if (x_tsc[i,j] == 0.0):
                    rhoenv[i,j] = 0.0
                if (x_tsc[i,j] > 0.0):
                    fv=(-1.0)*np.sqrt(xmo/x_tsc[i,j])
                    tsq=tau_tsc[i,j]*tau_tsc[i,j]
                    zeta=(xmo**3)*tsq/(16.*x_tsc[i,j])
                    ao=xmo/((x_tsc[i,j])**2.0)

                    # solve transcendental equation for theta_0
                    # first take care of special cases
                    # then, if theta_0 not yet assigned, call cubesolve
                    theta_0 = -999.0
                    pi2 = np.pi/2.

                    if zeta == 0.0:
                        theta_0 = tt[i,j]
                    if tt[i,j] == 0.0:
                        theta_0 = 0.0
                    if abs((tt[i,j]/np.pi)-1) <= 1.0e-5:
                        theta_0 = np.pi
                    if abs((tt[i,j]/pi2)-1) <= 1.0e-5:
                        if(zeta <= 1.0):
                            theta_0=pi2
                        else:
                            theta_0=-100.0

                    if theta_0 == -999.0: # continue on if needed
                        temp1=zeta
                        temp2=0.0
                        temp3=(1.0-zeta)
                        temp4=(-1.0)*np.cos(tt[i,j])
                        tsc_b=[temp1,temp2,temp3,temp4]
                        xone=1.0
                        result_cubesolve=cubesolve(tsc_b)

                        if result_cubesolve[3] > 0.0:
                            if (result_cubesolve[0] >= ((-1.0)*xone)) and (result_cubesolve[0] <= xone):
                                costheta_0=result_cubesolve[0]
                            else:
                                costheta_0=-0.5
                                print(f'Problem with cube root, costheta_0={result_cubesolve[0]}')

                        else:
                            costheta_0=-0.5
                            iroot=0
                            for z in [0,1,2]:
                                if ( (np.cos(tt[i,j])*result_cubesolve[z]) >= 0.0 and
                                      result_cubesolve[z] >= ((-1.0)*xone) and
                                      result_cubesolve[z] <= xone ):
                                    costheta_0=result_cubesolve[z]
                                    iroot=iroot+1

                        theta_0=np.arccos(costheta_0)

                    # now I have theta_0
                    ## begin handling special cases using the results of transcendental equation solution

                    # special case 1: in disk, on equator
                    rhoenv[i,j]=-999.0
                    if theta_0 < -1.0:
                        rhoenv[i,j]=0.0

                    # special case 2: outside disk, on equator
                    if abs((tt[i,j]/pi2)-1.0) <= 1.0e-5:
                        f4=1.+(2.*zeta*((1.-(1.5*np.sin(theta_0)))**(2.0)))
                        rhoenv[i,j]=((1./(4.*np.pi*gg*((modeltime*365.0*24.0*3600.0)**2.0)))
                                     *((-1.0)*ao)/(fv*sqrt(2.)*f4))

                    # general case if special cases not used
                    if rhoenv[i,j] == -999.0:
                        f1=np.sqrt(1.+(np.cos(tt[i,j])/np.cos(theta_0)))
                        f4=1.+(2.*zeta*((1.-(1.5*np.sin(theta_0)))**(2.0)))
                        rhoenv[i,j]=(1./(4.*np.pi*gg*((modeltime*365.0*24.0*3600.0)**2.0)))*((-1.0)*ao)/(fv*f1*f4)

                # switch off inner solution ## using...
                ## pre-set method for determining where the flattening becomes significant. outer solution = it hasn't
            if x_tsc[i,j] >= (tau_tsc[i,j]**2.0):  # TSC solution
                rhoenv[i,j]=((1./(4.*np.pi*gg*((modeltime*365.0*24.0*3600.0)**2.0)))
                             *(alpha_0_tsc[i,j]+((tau_tsc[i,j]**2.0)*(alpha_M_tsc[i,j]+(alpha_Q_tsc[i,j]*p2[i,j])))))
            # if rr[i,j] < renv_in:
            #     rhoenv[i,j]=rhofloor

            if not np.isfinite(rhoenv[i,j]):
                rhoenv[i,j]=0.0


    return rhoenv, u_r, y_tsc, x_tsc, tau_tsc[1,1]

def make_volumes(r_walls, theta_walls):
    rl  = r_integ(r_walls[:-1], r_walls[1:])            ## Integrate "length" between cell walls
    thl = th_integ(theta_walls[:-1], theta_walls[1:])   ## "
    rrl, ththl = np.meshgrid(rl, thl, indexing='ij') ## Mesh the "length" grids together
    cv = 2*np.pi*rrl*ththl  ## 'cell volume' of every cell in the simulation-space
    return cv

def trapezoid(x,f):
    if len(x) != len(f):
        return 'BAD FUNCTIONS'
    else:
        runsum = 0
        for i in range(1,len(nu)):
            dx = nu[i]-nu[i-1]
            py = fnu[i]+fnu[i-1]
            runsum += 0.5*dx*py
        return runsum

def bol_luminosity(nu,fnu): ## integrates a spectrum to get bolometric luminosity
    return trapezoid(nu, fnu)*4*np.pi*pc**2/ls

def bol_temperature(nu,fnu): ## calculates bolometric temperature (and confuses Colin with some weird quantum mechanics)
    return 1.25e-11*(trapezoid(nu, nu*fnu)/trapezoid(nu, fnu))

def bbody_lum(r, temp):
    return 4*np.pi*((r*rs)**2)*ssb*(temp**4)/ls

import numpy as np
from radmc_utils import *

cc = 29979245800 # speed of light cm/s

simdir = '/Users/cbourque/astrophysics/star-formation/SIMULATIONS'

get_modelnum = True
while get_modelnum:
  modelnumber = input('Generate input files for model number...')
  try:
    modelnumber = int(modelnumber)
    get_modelnum = False
  except:
    print('Input not an integer')
print(f'Calculating Lbol, Tbol for model{modelnumber:02n}')

mdir   = f'{simdir}/model{modelnumber:02n}'

mpars = np.loadtxt(f'{simdir}/model_parameters/model_parameters_{modelnumber}.tbl', skiprows=1) 
mpars = mpars[~np.isnan(mpars).any(axis=1)].T ## remove entries with nan values

times = mpars[0][1:]
lbols = []
tbols = []


for i in range(1,len(times)+1):
    intext   = np.loadtxt(f'{mdir}/OUTPUT/INTEXT/spectrum.{i:02n}.incl55.out', skiprows=2).T
    spec     = intext[1]
    nu       = 1e4*cc/intext[0]

    lb = bol_luminosity(nu, spec)
    tb = bol_temperature(nu, spec)

    lbols.append(lb)
    tbols.append(tb)

with open(f'lbol_tbol_model{modelnumber:02n}.tab', 'w+') as f:
    f.write(f'Model Time \t L_bol (L_sol) \t T_bol (K) \n')
    for i in range(1,len(times))):
        f.write(f'{1e6*times[i]:.6e} \t {lbols[i]:.6e} \t {tbols[i]:.6e} \n')

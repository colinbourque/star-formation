import numpy as np
from radmc_utils import *

cc = 29979245800 # speed of light cm/s

simdir = '/Users/cbourque/astrophysics/star-formation/SIMULATIONS'

get_modelnum = True
while get_modelnum:
  modelnumber = input('Calc indicators for model number...')
  try:
    modelnumber = int(modelnumber)
    get_modelnum = False
  except:
    print('Input not an integer')
print(f'Calculating Lbol, Tbol for model{modelnumber:02n}')

mdir   = f'{simdir}/model{modelnumber:02n}'

mpars = np.loadtxt(f'{simdir}/model_parameters/model_parameters_{modelnumber}.tbl', skiprows=1)
mpars = mpars[~np.isnan(mpars).any(axis=1)].T ## remove entries with nan values

times = mpars[0]
lbols = []
tbols = []


for i in range(1,len(mpars[0])):
    lb = []
    tb = []
    for inc in [5,15,25,35,45,55,65,75,85]:
        intext   = np.loadtxt(f'{mdir}/OUTPUT/INTEXT/spectrum.{i:03n}.incl{inc}.out', skiprows=2).T
        spec     = intext[1]
        nu       = 1e4*cc/intext[0]

        lb.append(bol_luminosity(nu, spec))
        tb.append(bol_temperature(nu, spec))

    lbols.append(lb)
    tbols.append(tb)

with open(f'lbols_model{modelnumber:02n}.tab', 'w+') as f:
    f.write(f'Model Time \t 5deg L_bol (L_sol) \t 15 \t 25 \t 35 \t 45 \t 55 \t 65 \t 75 \t 85 \n')
    for i in range(1,len(times)):
        f.write(f'{1e6*times[i]:.6e} \t {lbols[i][0]:.6e} \t {lbols[i][1]:.6e} \t {lbols[i][2]:.6e} \t {lbols[i][3]:.6e} \t {lbols[i][4]:.6e} \t {lbols[i][5]:.6e} \t {lbols[i][6]:.6e} \t {lbols[i][7]:.6e} \t {lbols[i][8]:.6e} \n')

with open(f'tbols_model{modelnumber:02n}.tab', 'w+') as f:
    f.write(f'Model Time \t 5deg T_bol (K) \t 15 \t 25 \t 35 \t 45 \t 55 \t 65 \t 75 \t 85 \n')
    for i in range(1,len(times)):
        f.write(f'{1e6*times[i]:.6e} \t {tbols[i][0]:.6e} \t {tbols[i][1]:.6e} \t {tbols[i][2]:.6e} \t {tbols[i][3]:.6e} \t {tbols[i][4]:.6e} \t {tbols[i][5]:.6e} \t {tbols[i][6]:.6e} \t {tbols[i][7]:.6e} \t {tbols[i][8]:.6e} \n')

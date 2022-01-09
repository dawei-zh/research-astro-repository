#!/usr/local/bin/python3

import numpy as np
import pandas as pd
from astropy.modeling import models,fitting
import healpy as hp

import warnings
warnings.filterwarnings("error")
errorlist = []

Norm = 10000
nside = 64
npixs = hp.nside2npix(nside)

v = models.Voigt1D()
fit = fitting.LevMarLSQFitter()

dfDark1 = pd.read_csv('IDark.csv')

# Step 1: Construct A Grid
part = 'I'
size = 0.1
rad = 45 # 45 for I, 300 for II, 600 for III
number = int(2*rad/size)
# Divide the whole space into different boxes
xbin = np.linspace(-rad, rad, number + 1)
ybin = np.linspace(-rad, rad, number + 1)
zbin = np.linspace(-rad, rad, number + 1)
# Put each particle into different bin
dfDark1['xbin'] = np.searchsorted(xbin, dfDark1['x'])
dfDark1['ybin'] = np.searchsorted(ybin, dfDark1['y'])
dfDark1['zbin'] = np.searchsorted(zbin, dfDark1['z'])

#Step 2: Use bin coordinate for LSR frame, so we could have velocity information for each bin
dfDark1['x'] = np.round(-rad + size*(dfDark1['xbin'] - 1) + size/2, decimals = 2) + 8.2
dfDark1['y'] = np.round(-rad + size*(dfDark1['ybin'] - 1) + size/2, decimals = 2)
dfDark1['z'] = np.round(-rad + size*(dfDark1['zbin'] - 1) + size/2, decimals = 2)
dfDark1['vx'] = dfDark1['vx'] - 20.380102
dfDark1['vy'] = dfDark1['vy'] - 224.709198
dfDark1['vz'] = dfDark1['vz'] - 3.895417

dfDark1['r'] = np.sqrt(dfDark1['x']**2+dfDark1['y']**2+dfDark1['z']**2)
dfDark1['vr'] = (dfDark1['x']*dfDark1['vx']+dfDark1['y']*dfDark1['vy']+dfDark1['z']*dfDark1['vz'])/dfDark1['r']

dfDark1.drop(['vx','vy','vz'], axis = 1, inplace = True)

#Step 3: Coordinate Transformation from galactic center frame to (l,b)
dfDark1['phi'] = np.arctan2(dfDark1['y'],dfDark1['x'])
dfDark1['theta'] = np.arccos(dfDark1['z']/dfDark1['r'])

dfDark1.drop(['r','x','y', 'z'], axis = 1, inplace = True)

dfDark1['hp'] = hp.ang2pix(nside, dfDark1['theta'], dfDark1['phi'])
dfDark1['E'] = -Norm*dfDark1['vr']/(300000+dfDark1['vr'])

print('I finished')

dfDark2 = pd.read_csv('IIDark.csv')

# Step 1: Construct A Grid

part = 'II'
size = 0.6
rad = 300 # 45 for I, 300 for II, 600 for III
number = int(2*rad/size)
# Divide the whole space into different boxes
xbin = np.linspace(-rad, rad, number + 1)
ybin = np.linspace(-rad, rad, number + 1)
zbin = np.linspace(-rad, rad, number + 1)
# Put each particle into different bin
dfDark2['xbin'] = np.searchsorted(xbin, dfDark2['x'])
dfDark2['ybin'] = np.searchsorted(ybin, dfDark2['y'])
dfDark2['zbin'] = np.searchsorted(zbin, dfDark2['z'])

#Step 2: Use bin coordinate for LSR frame, so we could have velocity information for each bin
dfDark2['x'] = np.round(-rad + size*(dfDark2['xbin'] - 1) + size/2, decimals = 2) + 8.2
dfDark2['y'] = np.round(-rad + size*(dfDark2['ybin'] - 1) + size/2, decimals = 2)
dfDark2['z'] = np.round(-rad + size*(dfDark2['zbin'] - 1) + size/2, decimals = 2)
dfDark2['vx'] = dfDark2['vx'] - 20.380102
dfDark2['vy'] = dfDark2['vy'] - 224.709198
dfDark2['vz'] = dfDark2['vz'] - 3.895417

dfDark2['r'] = np.sqrt(dfDark2['x']**2+dfDark2['y']**2+dfDark2['z']**2)
dfDark2['vr'] = (dfDark2['x']*dfDark2['vx']+dfDark2['y']*dfDark2['vy']+dfDark2['z']*dfDark2['vz'])/dfDark2['r']

dfDark2.drop(['vx','vy','vz'], axis = 1, inplace = True)

#Step 3: Coordinate Transformation from galactic center frame to (l,b)
dfDark2['phi'] = np.arctan2(dfDark2['y'],dfDark2['x'])
dfDark2['theta'] = np.arccos(dfDark2['z']/dfDark2['r'])

dfDark2.drop(['r','x','y', 'z'], axis = 1, inplace = True)

dfDark2['hp'] = hp.ang2pix(nside, dfDark2['theta'], dfDark2['phi'])
dfDark2['E'] = -Norm*dfDark2['vr']/(300000+dfDark2['vr'])

print('II finished')

res = pd.DataFrame(columns = ['hp','len','x0', 'A','fwhm_g', 'fwhm_l'])

steps = 0.8
for i in range(npixs):
    if i%2000 == 0:
        print('No. ' +str(i) + ' processing')
    Tmp1 = dfDark1[dfDark1['hp'] == i]
    Tmp2 = dfDark2[dfDark2['hp'] == i]
    length = len(Tmp1) + len(Tmp2)

    Tmp = Tmp1.append(Tmp2, ignore_index=True).copy()
    
    low = min(Tmp['E']) - 0.1
    high = max(Tmp['E']) + steps
    ebin = np.arange(start = low, stop = high, step = steps)
    
    Tmp['ebin'] = np.searchsorted(ebin, Tmp['E'])
    Tmp['ebin'] = low + steps/2 + (Tmp['ebin'] - 1) * steps
    Tmp.drop(['xbin','ybin','zbin', 'E', 'vr', 'phi', 'theta', 'hp'], axis=1, inplace=True)
    pdf = Tmp.groupby(by=['ebin']).sum().reset_index().copy()
    tot = np.sum(pdf['count'])
    pdf['count'] = pdf['count']/(tot*steps)

    
    try:
        tmpRes = fit(v,pdf['ebin'],pdf['count'],maxiter = 5000)
        res.loc[i] = [i, length, tmpRes.x_0.value, tmpRes.amplitude_L.value, tmpRes.fwhm_G.value, tmpRes.fwhm_L.value]
    except:
        errorlist.append("Voigt " + str(i) + "\n")
        Tmp.to_csv(str(i)+'.csv', index=False)
        res.loc[i] = [i, length, np.nan, np.nan, np.nan, np.nan]
        
    del Tmp
    del Tmp1
    del Tmp2
    
res.to_csv(r'./Voigt.csv', index=False)

wr = open('Voigt_err', 'a')
for i in range(len(errorlist)):
    wr.write(errorlist[i])

wr.close()




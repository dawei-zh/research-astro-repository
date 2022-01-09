#!/usr/local/bin/python3

import numpy as np
import pandas as pd

part = 'II'
size = 0.6
rad = 300 # 45 for I, 300 for II, 600 for III
number = int(2*rad/size)
# Divide the whole space into different boxes
xbin = np.linspace(-rad, rad, number + 1)
ybin = np.linspace(-rad, rad, number + 1)
zbin = np.linspace(-rad, rad, number + 1)
print('Step 1')
# Put each particle into different bin
dfDark = pd.read_csv('IIDark.csv')
dfDark['xbin'] = np.searchsorted(xbin, dfDark['x'])
dfDark['ybin'] = np.searchsorted(ybin, dfDark['y'])
dfDark['zbin'] = np.searchsorted(zbin, dfDark['z'])
print('Step 2')
# Calculate number of particles in each box
dfDark.drop(['x','y','z', 'vx', 'vy', 'vz'], axis = 1, inplace=True)
density = dfDark.groupby(['xbin','ybin','zbin']).sum()
density.reset_index(inplace=True)
print('Step 3')
# Estimate density for each box
density['rho'] = 35181.055*density['count']/(size**3)
density.drop(['count'], axis = 1, inplace = True)
print('Step 4')
# Central coordinate for each box under LSR Cartesian frame
density['x'] = np.round(-rad + size*(density['xbin'] - 1) + size/2, decimals = 2)
density['y'] = np.round(-rad + size*(density['ybin'] - 1) + size/2, decimals = 2)
density['z'] = np.round(-rad + size*(density['zbin'] - 1) + size/2, decimals = 2)
density.drop(['xbin','ybin','zbin'], axis = 1, inplace = True)
density.to_csv(r'./IIDensity.csv', index=False)
print('Step 5')

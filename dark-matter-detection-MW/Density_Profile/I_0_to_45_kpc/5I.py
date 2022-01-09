#!/usr/local/bin/python3

import numpy as np
import pandas as pd

part = 'I'
size = 0.1
rad = 45 # 45 for I, 300 for II, 600 for III
number = int(2*rad/size)
gridXYZ = pd.read_csv('gridIXYZ.csv')
print('Step 1')
density = pd.read_csv('IDensity.csv')
print('Step 2')
denGrid = pd.merge(gridXYZ, density, how='outer', on=['x','y', 'z'])
print('Step 3')
denGrid.fillna(0, inplace=True)
print('Step 4')
denGrid.to_csv(r'./IDenGrid.csv', index=False)
print('Step 5')

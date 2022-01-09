#!/usr/local/bin/python3

import numpy as np
import pandas as pd

part = 'I'
size = 0.1
rad = 45 # 45 for I, 300 for II, 600 for III
number = int(2*rad/size)
xGrid = np.round(np.linspace(-rad + size/2, rad - size/2, number), decimals = 2)
yGrid = np.round(np.linspace(-rad + size/2, rad - size/2, number), decimals = 2)
zGrid = np.round(np.linspace(-rad + size/2, rad - size/2, number), decimals = 2)
print('Step 1')
gridXYZ = pd.DataFrame(list(zGrid)*(number)*(number), columns = ['z'])
yy = np.array([])
xx = np.array([])
for i in range(number):
    tmpY = [yGrid[i]]*number
    yy = np.append(yy, tmpY)
    print(i)

yy = list(yy) * (number)
print('Step 2')
gridXYZ['y'] = yy
print('Step 3')
for i in range(number):
    tmpX = [xGrid[i]]*number*number
    xx = np.append(xx, tmpX)
    print(i)

gridXYZ['x'] = xx
print('Step 4')
gridXYZ.to_csv(r'./gridIXYZ.csv', index=False)
print('Step 5')

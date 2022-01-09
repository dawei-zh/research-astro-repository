#!/usr/local/bin/python3

# For new position np.pi/8, 3*np.pi/8, 5*np.pi/8, 7*np.pi/8, 9*np.pi/8, 11*np.pi/8, 13*np.pi/8, 15*np.pi/8
import numpy as np
import pandas as pd

part = 'II'
r_sc = 8.2 # 8.2 kpc
size = 0.6
rad = 300 # 45 for I, 300 for II, 600 for III
number = int(2*rad/size)
denGrid = pd.read_csv('IIDenGrid.csv')
print('Step 1')
xGrid = np.round(np.linspace(-rad + size/2, rad - size/2, number), decimals = 2)
yGrid = np.round(np.linspace(-rad + size/2, rad - size/2, number), decimals = 2)
zGrid = np.round(np.linspace(-rad + size/2, rad - size/2, number), decimals = 2)
Gridy, Gridx, Gridz = np.meshgrid(yGrid, xGrid, zGrid)
# Set an empty meshgrid for density 
rho = Gridx*0
flag = 0
for i in range(number):
    for j in range(number):
        for k in range(number):
            # rho[i][j][k] means density in different point
            # The iteration will extract density information from the top to the bottom
            # and the program technically works
            rho[i][j][k] = denGrid.loc[flag, 'rho']
            flag +=1
            if flag%10000 == 0:
                print(flag)

from scipy.interpolate import RegularGridInterpolator
profile = RegularGridInterpolator((xGrid, yGrid, zGrid), rho)
print('Step 2')

dftmp = pd.read_csv('hpLumin_128.csv')
dftmp.drop(['r3'], axis = 1, inplace = True)
error = []
lum = {}
lum2 = {}
Psi = np.array([np.pi/4, np.pi/2, 3*np.pi/4, np.pi, 5*np.pi/4, 3*np.pi/2, 7*np.pi/4])
for psi in Psi:
    lum2[int(360*psi/(2*np.pi))] = []
    lum[int(360*psi/(2*np.pi))] = []

def rot(xx, yy, zz, theta):
    x = xx*np.cos(theta) + yy*np.sin(theta)
    y = -xx*np.sin(theta) + yy*np.cos(theta)
    z = zz
    return [x,y,z]

print('Step 3')
for i in range(len(dftmp)):
    ll = dftmp.loc[i, 'long']
    bb = dftmp.loc[i, 'lat']
    ll = 2*np.pi*ll/360
    bb = 2*np.pi*bb/360
    r00 = dftmp.loc[i,'r1']//1
    r0 = dftmp.loc[i,'r2']

    r1 = (r0//0.001)/1000 - size/2 # round down to 3 decimals
    rnum = int((r1-r00)//size) + 2
    r = np.linspace(r00, r1, rnum)
    deltar = r[1] - r[0]
    for psi in Psi:
        rr = pd.DataFrame(r, columns=['r'])
        lsr = np.array([-r_sc*np.cos(psi), -r_sc*np.sin(psi)])
        ## 1. Galactic coordinate to cartesian
        rr['x'] = r * np.cos(bb) * np.cos(ll)
        rr['y'] = r * np.cos(bb) * np.sin(ll)
        rr['z'] = r * np.sin(bb)
        ## 2. Rotate counter clockwise psi
        rr['x'], rr['y'], rr['z'] = rot(rr['x'], rr['y'], rr['z'], -psi)
        ## 3. Transfer to Galactic center coordinate
        rr['x'] = rr['x'] + lsr[0]
        rr['y'] = rr['y'] + lsr[1]
        try:
            rr['lum2'] = profile(rr.loc[:,['x','y','z']])**2 * deltar
            rr['lum'] = profile(rr.loc[:,['x','y','z']]) * deltar
            lum2[int(360*psi/(2*np.pi))].append(np.log10(rr['lum2'].sum()))
            lum[int(360*psi/(2*np.pi))].append(np.log10(rr['lum'].sum()))
        except:
            error.append('Error '+str(i) + ' ' + str(int(360*psi/(2*np.pi))))
            lum2[int(360*psi/(2*np.pi))].append(np.nan)
            lum[int(360*psi/(2*np.pi))].append(np.nan)
        del rr
    if i%10000 == 0:
        print(i)



dftmp.drop(['r1', 'r2'], axis = 1, inplace = True)
for psi in Psi:
    dftmp['lum2 '+str(int(360*psi/(2*np.pi)))] = lum2[int(360*psi/(2*np.pi))]
    dftmp['lum '+str(int(360*psi/(2*np.pi)))] = lum[int(360*psi/(2*np.pi))]

dftmp.to_csv(r'./IIcurNew_128.csv', index=False)

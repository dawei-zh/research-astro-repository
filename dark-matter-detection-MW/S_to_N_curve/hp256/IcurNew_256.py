#!/usr/local/bin/python3

# For new position np.pi/8, 3*np.pi/8, 5*np.pi/8, 7*np.pi/8, 9*np.pi/8, 11*np.pi/8, 13*np.pi/8, 15*np.pi/8
import numpy as np
import pandas as pd
import healpy as hp

part = 'I'
r_sc = 8.2 # 8.2 kpc
size = 0.1
rad = 45 # 45 for I, 300 for II, 600 for III
number = int(2*rad/size)
denGrid = pd.read_csv('IDenGrid.csv')
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

dftmp = pd.read_csv('hpLumin.csv')
dftmp.drop(['r2', 'r3'], axis = 1, inplace = True)
error = []
lum = {}
lum2 = {}
Psi = np.array([np.pi/8, 3*np.pi/8, 5*np.pi/8, 7*np.pi/8, 9*np.pi/8, 11*np.pi/8, 13*np.pi/8, 15*np.pi/8])
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
    r0 = dftmp.loc[i,'r1']

    r1 = (r0//0.001)/1000 - 0.7*size # round down to 3 decimals
    rnum = int(r0//size) + 2
    r = np.delete(np.linspace(0, r1, rnum), 0)
    deltar = r[0]
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

dftmp.drop(['r1'], axis = 1, inplace = True)
for psi in Psi:
    dftmp['lum2 '+str(int(360*psi/(2*np.pi)))] = lum2[int(360*psi/(2*np.pi))]
    dftmp['lum '+str(int(360*psi/(2*np.pi)))] = lum[int(360*psi/(2*np.pi))]

dftmp.to_csv(r'./IcurNew_256.csv', index=False)

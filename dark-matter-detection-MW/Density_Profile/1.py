#!/usr/local/bin/python3

import numpy as np
import pandas as pd
import gizmo_read

# Load Simulation Data
dic = '.'
dmParticle = gizmo_read.read.Read.read_snapshot(species=['dark'], properties=['position', 'velocity'], directory=dic)

dmFilter = np.apply_along_axis(lambda x:np.sqrt(x[0]**2+x[1]**2+x[2]**2), 1, dmParticle['dark']['position'])

position = dmParticle['dark']['position']
velocity = dmParticle['dark']['velocity']
print('Step 0')

# Part I: Particles and their v <=45 kpc from GC
IDark = pd.DataFrame(position[dmFilter<=45], columns = ['x', 'y', 'z'])
vIDark = pd.DataFrame(velocity[dmFilter<=45], columns = ['vx', 'vy', 'vz'])
IDark['vx'] = vIDark['vx']
IDark['vy'] = vIDark['vy']
IDark['vz'] = vIDark['vz']
IDark['count'] = 1
IDark.to_csv(r'./IDark.csv', index=False)
del IDark
del vIDark
print('Step 1')

# Part II: Particles and their v 45<r<=300 kpc from GC
IIDark = pd.DataFrame(position[np.logical_and(dmFilter>45, dmFilter<=300)], columns = ['x', 'y', 'z'])
vIIDark = pd.DataFrame(velocity[np.logical_and(dmFilter>45, dmFilter<=300)], columns = ['vx', 'vy', 'vz'])
IIDark['vx'] = vIIDark['vx']
IIDark['vy'] = vIIDark['vy']
IIDark['vz'] = vIIDark['vz']
IIDark['count'] = 1
IIDark.to_csv(r'./IIDark.csv', index=False)
del IIDark 
del vIIDark
print('Step 2')

# Part III: Particles and their v 300<r<=600 kpc from GC
IIIDark = pd.DataFrame(position[np.logical_and(dmFilter>300, dmFilter<=600)], columns = ['x', 'y', 'z'])
vIIIDark = pd.DataFrame(velocity[np.logical_and(dmFilter>300, dmFilter<=600)], columns = ['vx', 'vy', 'vz'])
IIIDark['vx'] = vIIIDark['vx']
IIIDark['vy'] = vIIIDark['vy']
IIIDark['vz'] = vIIIDark['vz']
IIIDark['count'] = 1
IIIDark.to_csv(r'./IIIDark.csv', index=False)
print('Step 3')

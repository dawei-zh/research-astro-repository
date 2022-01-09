#!/usr/bin/python3.6

import numpy as np
import pandas as pd
from astropy.modeling import models,fitting
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("error")
errorlist = []

df = pd.read_csv('Voigt.csv')
npixs = 49152

steps = 0.8
xx = np.linspace(-30,30,60)
for i in range(npixs):
    if i%2000 == 0:
        print(i)
    para = df.iloc[i].values
    funVoigt = models.Voigt1D(para[2], para[3], para[5], para[4])    
    Tmp = pd.read_csv(r'/home/david/Documents/neutrino/Spectrum/6/csv/'+str(i)+'.csv')
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
        plt.figure()
        plt.plot(xx, funVoigt(xx))
        plt.scatter(pdf['ebin'], pdf['count'])
        plt.savefig(r'/home/david/Documents/neutrino/Spectrum/redo/pdf/'+str(i)+'.png')
        plt.close()
    except:
        errorlist.append("Voigt " + str(i) + "\n")





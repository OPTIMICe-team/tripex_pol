#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import matplotlib as mpl

import numpy as np
import xarray as xr

def maskData(variable, flag, mask):
    
    variable = variable.where((flag.values & mask) != mask)
    
    return variable

# dataPath = '/net/broebroe/lvonterz/tripexProcessing/output/20181201_tripex_pol_3fr_L2_mom.nc'

dataPath = '/work/lvonterz/tripexProcessing/output/20181101_tripex_pol_3fr_L2_mom.nc'
data = xr.open_dataset(dataPath)
print(data.range.values)
dbz_w = data['W_DBZ_H'].copy()
dbz_ka = data['Ka_DBZ_H'].copy()
qFlagW = data['quality_flag_offset_w'].copy()
print(data['offset_w'])
quit()
#masks 
maskP = int('1000000000000000',2) #n points
maskC = int('0100000000000000',2) # correl
maskV = int('0010000000000000',2) # Variance

#mask activation
# dbz_w = maskData(dbz_w, qFlagW, maskP)
dbz_w = maskData(dbz_w, qFlagW, maskC)
# dbz_w = maskData(dbz_w, qFlagW, maskV)
print(maskP)

plt.figure(figsize=(15,6))
(dbz_ka - dbz_w).plot(y='range', cmap='jet', vmin=-5, vmax=20)
plt.title('DWR_Ka_W')
plt.grid(b=True)
plt.show()



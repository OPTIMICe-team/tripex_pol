import glob
import numpy as np
import pandas as pd

import xarray as xr
from sys import argv
import glob
import os
from netCDF4 import Dataset
import time

inputPathW = '/data/obs/campaigns/tripex-pol/wband_gra/l0/'
date2proc = '20181208'
date = pd.to_datetime(date2proc)
filePathW = os.path.join(inputPathW,
                        date.strftime('%Y'),
                         date.strftime('%m'),
                         date.strftime('%d'))
filesW = glob.glob(filePathW+'/*LV0.nc')
Data = xr.open_dataset(filesW[0])

#Range = np.concatenate([Data.C1Range.values,Data.C2Range.values])

#for ind in range(len(Data.Time)):
#    Spec = np.concatenate([Data.C1VSpec[ind,:].values,Data.C2VSpec[ind,:].values])
#    if ind == 0:
#        VSpec = Spec
#    else:
#        VSpec = np.concatenate([VSpec,Spec])

C1Spec = xr.DataArray(Data.C1VSpec[1:2,:])
C1Spec = C1Spec.rename({'C1Range':'Range'})
C2Spec = xr.DataArray(Data.C2VSpec[1:2,:])
C2Spec = C2Spec.rename({'C2Range':'Range'})
C3Spec = xr.DataArray(Data.C3VSpec[1:2,:])
C3Spec = C3Spec.rename({'C3Range':'Range'})
C4Spec = xr.DataArray(Data.C4VSpec[1:2,:])
C4Spec = C4Spec.rename({'C4Range':'Range'})
Spec1 = C1Spec.combine_first(C2Spec)
#Spec2 = Spec1.combine_first(C3Spec)
#Spec = Spec2.combine_first(C4Spec)
print(Spec1)

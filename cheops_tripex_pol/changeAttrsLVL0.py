#-- this code changes the names and attributes from V to h polarization because that was asigned wrong before
import xarray as xr
import glob 
import pandas as pd

inputPath = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_0/'
inputFiles = sorted(glob.glob(inputPath+'*.nc'))
iDfile = '_tripex_pol_3fr_spec_filtered_regridded.nc'
outputPath = '/work/lvonterz/cheops_tripex_pol/output/LVL0/'
dateStart = pd.to_datetime('20181112'); dateEnd = pd.to_datetime('20190215')
dateList = pd.date_range(dateStart,dateEnd,freq='D')
for date in dateList:
    year = date.strftime('%Y')
    month = date.strftime('%m')
    day = date.strftime('%d')
    dateID = date.strftime('%Y%m%d')
    files = sorted(glob.glob(inputPath+year+'/'+month+'/'+day+'/*test.nc'))
    print(inputPath+year+'/'+month+'/spectra_x_filtered/'+day+'/*.nc')
    print(files)
    for f in files:
        print(f[-30:-28])
        hourID = f[-30:-28]
        data = xr.open_dataset(f)
        data = data.rename({'Ka_Z_V':'Ka_Z_H','Ka_VEL_V':'Ka_VEL_H',
                       'Ka_Z_H':'Ka_Z_V','Ka_VEL_H':'Ka_VEL_V',
                       'KaSpecV':'KaSpecH','KaSpecNoiseV':'KaSpecNoiseH',
                       'KaSpecH':'KaSpecV','KaSpecNoiseH':'KaSpecNoiseV'})
        data = data.rename({'W_Z_V':'W_Z_H','W_VEL_V':'W_VEL_H',
                       'W_Z_H':'W_Z_V','W_VEL_H':'W_VEL_V',
                       'WSpecV':'WSpecH','WSpecNoiseV':'WSpecNoiseH',
                       'WSpecH':'WSpecV','WSpecNoiseH':'WSpecNoiseV'})
        data = data.rename({'X_Z_V':'X_Z_H','X_VEL_V':'X_VEL_H',
                        'XSpecV':'XSpecH','XSpecNoiseV':'XSpecNoiseH'})

        data.Ka_Z_H.attrs['long_name'] = 'equivalent reflectivity factor, horizontal polarization, Ka Band'
        data.Ka_Z_V.attrs['long_name'] = 'equivalent reflectivity factor, vertical polarization, Ka Band'
        data.W_Z_H.attrs['long_name'] = 'equivalent reflectivity factor, horizontal polarization, W Band'
        data.W_Z_V.attrs['long_name'] = 'equivalent reflectivity factor, vertical polarization, W Band'
        data.X_Z_H.attrs['long_name'] = 'equivalent reflectivity factor, horizontal polarization, X Band'

        data.Ka_VEL_H.attrs['long_name'] = 'mean Doppler velocity, horizontal polarization, Ka Band'
        data.Ka_VEL_V.attrs['long_name'] = 'mean Doppler velocity, vertical polarization, Ka Band'
        data.W_VEL_H.attrs['long_name'] = 'mean Doppler velocity, horizontal polarization, W Band'
        data.W_VEL_V.attrs['long_name'] = 'mean Doppler velocity, vertical polarization, W Band'
        data.X_VEL_H.attrs['long_name'] = 'mean Doppler velocity, horizontal polarization, X Band'

        data.KaSpecH.attrs['long_name'] = 'Doppler Spectrum, horizontal polarization, Ka Band'
        data.KaSpecV.attrs['long_name'] = 'Doppler Spectrum, vertical polarization, Ka Band'
        data.WSpecH.attrs['long_name'] = 'Doppler Spectrum, horizontal polarization, W Band'
        data.WSpecV.attrs['long_name'] = 'Doppler Spectrum, vertical polarization, W Band'
        data.XSpecH.attrs['long_name'] = 'Doppler Spectrum, horizontal polarization, X Band'

        data.KaSpecNoiseH.attrs['long_name'] = 'spectral Noise Level, horizontal polarization, Ka Band'
        data.KaSpecNoiseV.attrs['long_name'] = 'spectral Noise Level, vertical polarization, Ka Band'
        data.WSpecNoiseH.attrs['long_name'] = 'spectral Noise Level, horizontal polarization, W Band'
        data.WSpecNoiseV.attrs['long_name'] = 'spectral Noise Level, vertical polarization, W Band'
        data.XSpecNoiseH.attrs['long_name'] = 'spectral Noise Level, horizontal polarization, X Band'

        data.to_netcdf(outputPath+'/'+year+'/'+month+'/'+day+'/'+dateID+'_'+hourID+iDfile)

    quit()




















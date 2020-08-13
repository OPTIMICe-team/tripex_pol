#-- this code changes the names and attributes from V to h polarization because that was asigned wrong before
import xarray as xr
import glob 

inputPath = '/work/lvonterz/cheops_tripex_pol/output/'
inputFiles = sorted(glob.glob(inputPath+'*moments.nc'))
iDfile = '_tripex_pol_3fr_L1_mom.nc'
print(inputFiles)
for f in inputFiles:
    print(f)
    print(f[-19:-11])
    date = f[-19:-11]
    print(date)
    data = xr.open_dataset(f)
    data = data.rename({'Ka_DBZ_V':'Ka_DBZ_H','Ka_VEL_V':'Ka_VEL_H',
                       'Ka_WIDTH_V':'Ka_WIDTH_H','Ka_SK_V':'Ka_SK_H',
                       'Ka_DBZ_H':'Ka_DBZ_V','Ka_VEL_H':'Ka_VEL_V',
                       'Ka_WIDTH_H':'Ka_WIDTH_V','Ka_SK_H':'Ka_SK_V'})
    data = data.rename({'W_DBZ_V':'W_DBZ_H','W_VEL_V':'W_VEL_H',
                       'W_WIDTH_V':'W_WIDTH_H','W_SK_V':'W_SK_H',
                       'W_DBZ_H':'W_DBZ_V','W_VEL_H':'W_VEL_V',
                       'W_WIDTH_H':'W_WIDTH_V','W_SK_H':'W_SK_V'})
    data = data.rename({'X_DBZ_V':'X_DBZ_H','X_VEL_V':'X_VEL_H',
                       'X_WIDTH_V':'X_WIDTH_H','X_SK_V':'X_SK_H'})

    data.Ka_DBZ_H.attrs['long_name'] = 'equivalent reflectivity factor, horizontal polarization, Ka Band'
    data.Ka_DBZ_V.attrs['long_name'] = 'equivalent reflectivity factor, vertical polarization, Ka Band'
    data.W_DBZ_H.attrs['long_name'] = 'equivalent reflectivity factor, horizontal polarization, W Band'
    data.W_DBZ_V.attrs['long_name'] = 'equivalent reflectivity factor, vertical polarization, W Band'
    data.X_DBZ_H.attrs['long_name'] = 'equivalent reflectivity factor, horizontal polarization, X Band'

    data.Ka_VEL_H.attrs['long_name'] = 'mean Doppler velocity, horizontal polarization, Ka Band'
    data.Ka_VEL_V.attrs['long_name'] = 'mean Doppler velocity, vertical polarization, Ka Band'
    data.W_VEL_H.attrs['long_name'] = 'mean Doppler velocity, horizontal polarization, W Band'
    data.W_VEL_V.attrs['long_name'] = 'mean Doppler velocity, vertical polarization, W Band'
    data.X_VEL_H.attrs['long_name'] = 'mean Doppler velocity, horizontal polarization, X Band'

    data.Ka_WIDTH_H.attrs['long_name'] = 'radar spectral width, horizontal polarization, Ka Band'
    data.Ka_WIDTH_V.attrs['long_name'] = 'radar spectral width, vertical polarization, Ka Band'
    data.W_WIDTH_H.attrs['long_name'] = 'radar spectral width, horizontal polarization, W Band'
    data.W_WIDTH_V.attrs['long_name'] = 'radar spectral width, vertical polarization, W Band'
    data.X_WIDTH_H.attrs['long_name'] = 'radar spectral width, horizontal polarization, X Band'

    data.Ka_SK_H.attrs['long_name'] = 'skewness, horizontal polarization, Ka Band'
    data.Ka_SK_V.attrs['long_name'] = 'skewness, vertical polarization, Ka Band'
    data.W_SK_H.attrs['long_name'] = 'skewness, horizontal polarization, W Band'
    data.W_SK_V.attrs['long_name'] = 'skewness, vertical polarization, W Band'
    data.X_SK_H.attrs['long_name'] = 'skewness, horizontal polarization, X Band'

    data.to_netcdf('../output/'+date+iDfile)






















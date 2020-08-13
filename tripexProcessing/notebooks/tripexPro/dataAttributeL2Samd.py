def globalAttributes(rootgrpOut):
   
   rootgrpOut.Title = 'Regridded and corrected ground-based X, Ka, W-Band Doppler cloud radar moments obtained during the TRIple-frequency and Polarimetric radar Experiment for improving process observation of winter precipitation (TRIPEx) campaign at JOYCE supersite (Research Center Juelich)'
   rootgrpOut.Institution = 'Emmy-Noether Group OPTIMIce, Institute for Geophysics and Meteorology (IGMK), University of Cologne'
   rootgrpOut.Contact_person = 'Leonie von Terzi, lterzi@uni-koeln.de'
   rootgrpOut.Source = 'METEOR-50-DX (X-Band radar owned by KIT), Metek MIRA-35 (Ka-Band radar JOYRAD-35 owned by IGMK), RPG-FMCW-94-SP (W-Band radar JOYRAD-94 owned by IGMK)' 
   rootgrpOut.History = 'Data processed with tripexPro.py by University of Cologne'
#   rootgrpOut.Dependencies = 'External'
   rootgrpOut.Conventions = 'CF-1.6 where applicable'
   rootgrpOut.Author = 'Jose Dias Neto (jdiasnet@uni-koeln.de)'
   rootgrpOut.Comments = 'The original radar moment data were re-gridded in height and time; attenuation and offset corrections have been applied but can be reversed by the user since the corrections applied are stored; polarimetry data from  Ka Band (LDR) and from X Band (PhiDP, ZDR, RhoHV) is included. No additional processing was applied to PhiDP, ZDR and RhoHV (They are output from the radar software).'
#   rootgrpOut.License = 'For non-commercial use only. This data is subject to the HD(CP)2 data policy to be found at hdcp2.zmaw.de and in the HD(CP)2 Observation Data Product Standard'
   rootgrpOut.Pulse_repetition_frequency_x = '1.2 kHz'
   rootgrpOut.Number_of_FFT_x = '1200'
   rootgrpOut.Number_of_spectral_average_x = '1'
   rootgrpOut.Pulse_repetition_frequency_ka = '5 kHz'
   rootgrpOut.Number_of_FFT_ka = '512'
   rootgrpOut.Number_of_spectral_average_ka = '20'
   rootgrpOut.Chirp_squence_w = 'Chirp sequence: 1st (from 100 to 400 m); 2nd (from 400 to 1200 m);  3rd (from 1200 to 3000 m);  4th (from 3000 to 12000 m)'
   rootgrpOut.Pulse_repetition_frequency_w = 'For each chirp sequence: 1st (12.2 kHz); 2nd (10.2 kHz); 3rd (7.8 kHz); 4th (5.3 kHz)'
   rootgrpOut.Number_of_FFT_w = 'For each chirp sequence: 1st (512); 2nd (512); 3rd (512); 4th (512)'
   rootgrpOut.Number_of_spectral_average_w = 'For each chirp sequence: 1st (8); 2nd (8); 3rd (8); 4th (18)'
   return rootgrpOut


def timeAttributes(timeRef):
  
   timeRef.standard_name = 'time'
   timeRef.long_name = 'time'
   timeRef.comments = 'center of averaging period'  
   timeRef.units = 'seconds since 1970-01-01 00:00:00 UTC'
   timeRef.bounds = 'time_bnds'

   return timeRef


def rangeAttributes(rangeRef):

   rangeRef.long_name = 'distance from sensor to center of each range gates along the line of sight'
   rangeRef.comments = 'Vertical distance (because zenith operation only) of the center radar range gates to the antenna of the reference radar (W band)'
   rangeRef.units = 'm'
   rangeRef.bounds = "range_bnds" 

   return rangeRef


def variableAttribute(variable, varName, radar):
  
   radar = radar.capitalize()
   if radar == 'X':
      commentZe = 'The original reflectivities have been re-gridded onto a common time-height grid and have been corrected for gas attenuation; relative calibration has been applied using Ka-Band radar as reference in the high altitude ice part when possible. Note that the reflectivity resolution of the instrument is only 0.5 dB.'

      commentVel = 'Sign convention: Negative when moving towards the radar. Note that the radar antenna had to continuously rotate (bird bath scan) in order to record data. The antenna seemed to had a slight mis-pointing thus a sinusoidal modulation is sometimes visible depending on the horizontal wind. This signature disappears when applying a temporal averaging of 3 minutes which matches the duration for a complete antenna rotation.'

      commentOff = 'After gas attenuation correction, the remaining reflectivity offsets due to wet radome, liquid water attenuation, etc. have been estimated by using a 15 min window, regions at least 1km above the zero degree isotherm and reflectivities smaller than -10dBZ. Assuming the Ka-band radar as our reference, we assume that for X and Ka-band these low-reflectivity ice clouds represent Rayleigh scatterers and hence the reflectivities should match. A quality flag is provided indicating the reliability of the correction applied. The user can simply go back to the attenuation-only corrected reflectivities by subtracting this offset from the corrected reflectivity values'

   elif radar == 'Ka':
      commentZe = 'The original reflectivities have been re-gridded onto a common time-height grid and have been corrected for gas attenuation; relative calibration has been performed with rain cases observed together with PARSIVEL.'
      
      commentVel = 'Sign convention: Negative when moving towards the radar. Note that we are suspecting a slight mispointing of the antenna in the order of 0.5 degree off-zenith. This needs more in-depths investigation but can lead to offsets in the order of 0.1-0.3 m/s in case of strong horizontal winds.'

   elif radar == 'W':
      commentZe = 'The original reflectivities have been re-gridded onto a common time-height grid and a Ghost-Echo filtering has been applied. After, the reflectivities have been corrected for gas attenuation; relative calibration has been applied using Ka-Band radar as reference in the high altitude ice part when possible.'
      
      commentVel = 'Sign convention: Negative when moving towards the radar'

      commentOff = 'After gas attenuation correction, the remaining reflectivity offsets due to wet radome, liquid water attenuation, etc. have been estimated by using a 15 min window, regions at least 1km above the zero degree isotherm and reflectivities smaller than -10dBZ. Assuming the Ka-band radar as our reference, we assume that for W and Ka-band these low-reflectivity ice clouds represent Rayleigh scatterers and hence the reflectivities should match. A quality flag is provided indicating the reliability of the correction applied. The user can simply go back to the attenuation-only corrected reflectivities by subtracting this offset from the corrected reflectivity values'

   if varName == 'DBZ':
      standard_name = 'equivalent_reflectivity_factor'
      long_name = 'equivalent reflectivity factor {radar}-Band'.format(radar=radar)
      units = 'dBZ'
      finalComment = commentZe
   
   elif varName == 'VEL':
      standard_name = 'radial_velocity_of_scatterers_away_from_instrument'
      long_name = 'mean Doppler velocity {radar}-Band'.format(radar=radar)
      units = 'm s-1'
      finalComment = commentVel

   elif varName == 'WIDTH':
      long_name = 'radar spectral width {radar}-Band'.format(radar=radar)
      units = 'm s-1'

   elif varName == 'SK':
      long_name = 'radar skewness {radar}-Band'.format(radar=radar)
 
   elif varName == 'ldr':
      long_name = 'radar linear depolarization ratio {radar}-Band'.format(radar=radar)
      units = 'dB'

   elif varName == 'kdp':
      long_name = 'radar specific differential phase {radar}-Band'.format(radar=radar)
      units = 'deg km-1'

   elif varName == 'phidp':
      long_name = 'radar integrated differential phase {radar}-Band'.format(radar=radar)
      units = 'deg'
 
   elif varName == 'zdr':
      long_name = 'radar differential reflectivity {radar}-Band'.format(radar=radar)
      units = 'dB'

   elif varName == 'rhohv':
      long_name = 'radar co-polar correlation function {radar}-Band'.format(radar=radar)
      units = '1'
  
   elif varName == 'pia':
      long_name = 'path integrated attenuation {radar}-Band (two-way) due to gaseous atmosphere'.format(radar=radar)
      units = 'dB'
      finalComment = 'Two-way path integrated attenuation due to dry gases and water vapor calculated with radiative transfer model PAMTRA and thermodynamic profiles from CloudNet product for JOYCE site. Simply subtract these values from the corrected reflectivities in order to restore the uncorrected reflectivities.'

   elif varName == 'offset':
      long_name = 'relative offset of {radar}-Band reflectivities compared to reference Ka-Band'.format(radar=radar)
      units = 'dBZ'
      finalComment = commentOff

   elif varName == 'freq_sb':
      standard_name = 'sensor_band_central_radiation_frequency'
      long_name = 'operating frequency of {radar}-Band radar'.format(radar=radar)
      units = 's-1'

   elif varName == 'radar_beam_width':
      long_name = '3dB beam width of the {radar}-Band radar'.format(radar=radar)
      units = 'degree'

   elif varName == 'lon':
      standard_name = 'longitude'	
      long_name = 'longitude of reference radar system' 
      finalComment = 'Longitude of the reference radar (W-Band) location in decimals.'
      units = 'degrees_east'

   elif varName == 'lat':
      standard_name = 'latitude'	
      long_name = 'latitude of reference radar system' 
      finalComment = 'Latitude of the reference radar (W-Band) location in decimals.'
      units = 'degrees_north'

   elif varName == 'ta':
      standard_name = 'air_temperature'	
      long_name = 'temperature' 
      source = 'cloudnet model'
      units = 'deg C'

   elif varName == 'pa':
      standard_name = 'air_pressure'	
      long_name = 'pressure' 
      source = 'cloudnet model'
      units = 'Pa'

   elif varName == 'hur':
      standard_name = 'relative_humidity'	
      long_name = 'relative humidity' 
      source = 'cloudnet model'
      units = '1'


   elif varName == 'zsl':
      standard_name = 'altitude' 
      long_name = 'altitude above mean sea level'
      finalComment = 'altitude above mean sea level of the reference radar (W-Band) location'
      units = 'm'

   elif varName == 'quality_flag_offset':
      long_name = 'quality flag indicating reliability of the offset correction for the {radar}-Band data'.format(radar=radar)
      finalComment = 'Bits 0 to 5: empty; Bit 6: 0 if the liquid water path is less than 200 g m-2, 1 if the liquid water path is greater than 200 g m-2; Bit 7: 0 no rain, 1 rain; Bits 8 to 12: empty; Bit 13: 0 if the variance of the DWR_{radar}_Ka within a 3 minutes time window is less than 2 dB**2, 1 if the DWR variance is greater than 2 dB**2; Bit 14: 0 if correlation of {radar} and Ka Band reflectivities is larger than 0.7, 1 if correlation is less than 0.7; Bit 15: 0 if the number of points used to calculate the offset is greater than 300, 1 if the number of points is less than 300. If Bit 14 or higher is set, we recommend not to use the calculated offsets but e.g. rather interpolate between time periods with high-quality offset estimates.'.format(radar=radar)
      units = '1'

   try: 
        variable.standard_name = standard_name
   except:
        print('no standard_name to write')
        pass
  
   try:
        variable.long_name = long_name
   except:
        print('no long_name to write')
        pass

   try:
        variable.units = units
   except:
        print('no units to write')
        pass

   try:
        variable.comments = finalComment
   except:
        print('no comment to write')
        pass

   return variable





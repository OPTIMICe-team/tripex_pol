def globalAttributes(rootgrpOut):
   
   rootgrpOut.Experiment = 'TRIPEX, Forschungszentrum Juelich'
   rootgrpOut.Instruments = 'JOYRAD-94 94.4 GHz cloud radar (W band), JOYRAD-35 35.5 GHz cloud radar (Ka band), KIXPOL 9.4 GHz radar (KIT) (X band)'
   rootgrpOut.Data = 'Produced by Jose Dias, jdiasnet@uni-koeln.de'
   rootgrpOut.Processed_with = 'tripexPro*.py'
   rootgrpOut.Institution = 'Data processed within the Emmy-Noether Group OPTIMIce, Institute for Geophysics and Meteorology, University of Cologne, Germany'
   rootgrpOut.comment = 'All radar data have been only resampled for this data product. No offset or attenuation corrections applied'
   rootgrpOut.Latitude = '50.908547 N'
   rootgrpOut.Longitude = '6.413536 E'
   rootgrpOut.Altitude = 'Altitude of the JOYCE (www.joyce.cloud) platform: 111m asl'

   return rootgrpOut


def timeAttributes(timeRef):
   
   timeRef.long_name = 'time in sec since 01.01.1970 00:00:00'
   timeRef.units = 's'

   return timeRef


def rangeAttributes(rangeRef):

   rangeRef.long_name = 'Vertical distance of the center radar range gates to the JOYCE platform + 1.5m (altitude of W band center range gates)'
   rangeRef.units = 'm'

   return rangeRef


def variableAttribute(variable, varName, radar):

   if radar == 'X':
      source = 'KIXPOL 9.4 GHz radar (KIT)'
      velocityLim = '+/- 80 m s-1'

   if radar == 'Ka':
      source = 'JOYRAD-35 35.5 GHz cloud radar'
      velocityLim = '+/- 20 m s-1'

   if radar == 'W':
      source = 'JOYRAD-94 94 GHz cloud radar'
      velocityLim = '+/- 18 m s-1'

   if varName == 'Ze':
      long_name = 'Equivalent '+radar+' band Reflectivity Factor Ze of all Targets'
      units = 'dBZ'

   if varName == 'Ze' and radar == 'Ka':
      offset_correction = 'The Joyrad-35 reflectivities obtained from the mmclx files have to be corrected: 1) An offset of 2dB has to be added in order to compensate for the finite receiver beamwidth loss. 2) A recently discovered offset of additional 3dB has to be added to compensate for some errors in the DSP processing. For further information, contact Metek (bauer@metek.de).'
      variable.offset_correction = offset_correction     

   if varName == 'v':
      long_name = radar+' band Mean Doppler velocity (Sign convention: Negative when moving towards the radar) '
      units = 'm s-1'
      variable.Nyquist_velocity = velocityLim

   if varName == 'SW':
      long_name = radar+' band Spectrum Width'
      units = 'm s-1'
 
   if varName == 'LDR':
      long_name = radar+' band Linear De-Polarization Ratio'
      units = 'dB'
    
   if varName == 'delta_altitude':
      long_name = radar+' band vertical distance of the original range resolution to the vertical grid (var:altitude) used to resample the data'
      units = 'm'

   if varName == 'delta_time':
      long_name = radar+' band temporal difference between the original time resolution and the time vector used to resample the data'
      units = 's'
 
     
   variable.long_name = long_name
   variable.units = units
   variable.source = source

   return variable





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 01:03:00 2019

@author: Delano Yoder
"""
# Timer
import time
start_timer = time.time()

import math
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# File name
file_name = 'LMC_ATCA+PKS_K_integrated_trimmed_rebinned.fits'

# Importing fits header and data
fits_file = fits.open(file_name)
header, data = fits_file[0].header, fits_file[0].data

# Integrating and normalizing data into an image
while (len(data.shape) > 2):
    data = np.nansum(data, 0)
image = data/np.amax(data[np.where(~np.isnan(data))])

# Variables for image width and length
image_width, image_length = image.shape[0], image.shape[1]

# Determining half the length of the shortest side of the image
max_r = image_length/2.0 if image_length <= image_width else image_width/2.0


# Structure function dictionary variable
SF_dict = {}

# Setting variables for structure function
w, l = image_width, image_length

# Switch for percentage print out
switch = True

# Stepping through all different combinations of pixel separations
for x in range(0, int(max_r)):
    for y in range(0, int(max_r)):
        
        # Distance between the two arbitrary pixels
        r = (x**2+y**2)**0.5
        
        # If the pixel separation does not exceed the maximum limit
        if (r <= max_r):
            
            # Calculate all pixel separations
            SFx = np.nanmean((image[0:w-x,0:l-y]-image[x:w,y:l])**2)
            SFy = np.nanmean((image[x:w,0:l-y]-image[0:w-x,y:l])**2)
        
            # Storing structure function values based on pixel separation
        
            # Avoiding double calculations parallel movements
            if (x == 0 or y == 0):
                if (r in SF_dict):
                    SF_dict[r].append(SFx)
                else:
                    SF_dict[r] = [SFx]
                
            # No double calculations on diagonal movements
            else:
                if (r in SF_dict):
                    SF_dict[r].append(SFx)
                    SF_dict[r].append(SFy)
                else:
                    SF_dict[r] = [SFx, SFy]
                    
        # Time estimate
        if (x*y > (max_r/2)**2 and switch != False):
            time_estimate = (time.time()-start_timer)*2.86
            if (time_estimate  <= 60):
                print('Estimated time is '+str(time_estimate)+'secs')
            elif (time_estimate > 60 and time_estimate <= 3600):
                print('Estimated time is '+str(time_estimate/60.)+'mins')
            else:
                print('Estimated time is '+str(time_estimate/3600.)+'hrs')
            print('25%')
            switch = False
            
                    
# Removing the pixel separation of 0
SF_dict.pop(0)
            
# Allocating memory for pixel separation array, structure function average
# array, and standard deviation array
lags = np.empty(len(SF_dict.keys()))*np.nan
avgs = np.empty(len(SF_dict.keys()))*np.nan
stds = np.empty(len(SF_dict.keys()))*np.nan

# Calculating the mean and standard deviation of each structure function 
# respectively based on pixel separation
index = 0
for key,value in SF_dict.items():
    lags[index] = key
    avgs[index] = np.nanmean(value)
    stds[index] = np.nanstd(value)
    index += 1
    
# Sorting structure function bins of spaced by 0.05 pixels logarithmically    
x_SF = []
y_SF = []
e_SF = []
numBins = math.ceil(np.log10(max_r)/0.05)
for bins in range(numBins):
    ind = (np.log10(lags) >= bins*0.05)*(np.log10(lags) <= (bins+1)*0.05)
    if (True in ind):
        x_SF.append(np.nanmean(lags[ind]))
        y_SF.append(np.nansum(avgs[ind]/stds[ind]**2)/np.nansum(1/stds[ind]**2))
        e_SF.append((1/np.nansum(1/stds[ind]**2))**0.5)
        
# Formatting structure function results into a list
SF = []       
SF.append(np.asarray(x_SF))
SF.append(np.asarray(y_SF))
SF.append(np.asarray(e_SF))

# Variable for new fits file name
try:
    fits_name = header['OBJECT'][:3].upper()+'_'+header['TELESCOP'].upper()+'_INT+SF'
except KeyError:
    fits_name = file_name[:-5]+'_INT+SF'
    
# Creating the new fits file containing the structure function results
try:
    fits.PrimaryHDU(SF, header).writeto(fits_name+'.fits')
except IOError:
    rn = '%.4f' % np.random.rand(1)[0]
    print('Fits files arleady exist.')
    print('New file name is '+rn)
    fits.PrimaryHDU(SF, header).writeto(fits_name+rn+'.fits')
    
# Plotting the structure function results
plt.figure(figsize=(5, 5))
plt.axes(xscale='log')
plt.xlabel('Linear Scale [pc]')
plt.ylabel("Log(SF(r'))")
plt.plot(SF[0], np.log10(SF[1]), 'bs', markersize=7, markerfacecolor='None')
plt.show()
        
print(str(time.time()-start_timer)+'s')
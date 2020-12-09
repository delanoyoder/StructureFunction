#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 12:11:11 2019

@author: Delano Yoder
"""
# Timer
import time
start_timer = time.time()

import math
import numpy as np
from astropy.io import fits

# File name
file_name = 'smc_atca+pks.imbin.fits'

# Importing fits header and data
fits_file = fits.open(file_name)
header, data = fits_file[0].header, fits_file[0].data

# Check if the data has enough dimensions
if (len(data.shape) < 3):
    print("Data does not have enough dimensions")

# Integrating and normalizing data into an image
while (len(data.shape) > 3):
    data = np.nansum(data, 0)
image = data/np.amax(data[np.where(~np.isnan(data))])

# Variables for image width, length and depth
image_depth, image_width, image_length = image.shape[0], image.shape[1], image.shape[2]

# Kernel diameter and step size
kernel_diameter, step_size = 51, 1

# Number of possible kernel steps across the image and down the image
number_x_kernels = int((image_width-kernel_diameter)/step_size+1)
number_y_kernels = int((image_length-kernel_diameter)/step_size+1)

# Structure Function
def SF(image):
    
    # Storing the length and width of each image into variables
    image_length, image_width = image.shape[0], image.shape[1]
    
    # Determining half the length of the shortest side of the image
    max_r = image_length/2.0 if image_length <= image_width else image_width/2.0

    # Structure function dictionary variable
    SF_dict = {}

    # Setting variables for structure function
    l, w = image_length, image_width

    # Stepping through all different combinations of pixel separations
    for x in range(0, int(max_r)):
        for y in range(0, int(max_r)):
        
            # Distance between the two arbitrary pixels
            r = (x**2+y**2)**0.5
        
            # If the pixel separation does not exceed the maximum limit
            if (r <= max_r):
            
                # Calculate all pixel separations
                SFx = np.nanmean((image[0:l-x,0:w-y]-image[x:l,y:w])**2)
                SFy = np.nanmean((image[x:l,0:w-y]-image[0:l-x,y:w])**2)
        
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
            
    return np.asarray(x_SF), np.asarray(y_SF), np.asarray(e_SF)

def RSF(image):
    # Stepping through each kernel
    for x in range(number_x_kernels):
        for y in range(number_y_kernels):
            
            # Variables for the position of each edge of the kernel
            top = y*step_size
            left = x*step_size
            bottom = y*step_size+kernel_diameter
            right = x*step_size+kernel_diameter
            
            # Checking if kernel is within the image
            if (right <= image_width and bottom <= image_length):
                kernel = image[left:right, top:bottom].copy()
                
                # Calculating maximum limit
                max_r = kernel_diameter/2.0
                
                # Making the kernal circlular   
                r_2 = (np.arange(0, kernel_diameter, dtype=float)-max_r)**2
                cdist = np.empty((kernel_diameter, kernel_diameter))
                for i in range(kernel_diameter):
                    cdist[0:, i] = (r_2+r_2[i])**0.5
                kernel[cdist > max_r] = np.nan
                
                # Performing structure function on the kernel
                x_SF, y_SF, e_SF = SF(kernel)
                
                # Allocating memory for the results
                if (x+y == 0):
                    Y = np.empty((number_x_kernels, number_y_kernels, len(y_SF)))*np.nan
                    E = np.empty((number_x_kernels, number_y_kernels, len(e_SF)))*np.nan
                
                # Storing the results
                Y[x, y, :] = y_SF
                E[x, y, :] = e_SF
                
    return x_SF, Y, E

# Stepping through each channel
for c in range(image_depth):
    
    # Applying the structure function to an individual channel
    X, y, e = RSF(image[c,:,:])

    # Allocating memory for the structure function in each channel
    if (c == 0):
        Y = np.empty((image_depth, y.shape[0], y.shape[1], y.shape[2]))*np.nan
        E = np.empty((image_depth, e.shape[0], e.shape[1], e.shape[2]))*np.nan
        
        # Time estimate
        time_estimate = (time.time()-start_timer)*image_depth
        if (time_estimate  <= 60):
            print('Estimated time is '+str(time_estimate)+'secs')
        elif (time_estimate > 60 and time_estimate <= 3600):
            print('Estimated time is '+str(time_estimate/60.)+'mins')
        else:
            print('Estimated time is '+str(time_estimate/3600.)+'hrs')
            
    if (c%8 == 0):   
        print(str(100.*((c+1)/image_depth))+'%')

    # Storing the channel's structure function results
    Y[c, :, :, :] = y
    E[c, :, :, :] = e

# Variable for new fits file name
try:
    fits_name = header['OBJECT'][:3].upper()+'_'+header['TELESCOP'].upper()
    fits_name += '_SCA+RSF_'+str(kernel_diameter)+'_'+str(step_size)
except KeyError:
    fits_name = file_name[:-5]+'_SCA+RSF_'+str(kernel_diameter)+'_'+str(step_size)
    
# Creating the new fits file containing the structure function results
try:
    fits.PrimaryHDU(X, header).writeto(fits_name+'_X.fits')
    fits.PrimaryHDU(Y, header).writeto(fits_name+'_Y.fits')
    fits.PrimaryHDU(E, header).writeto(fits_name+'_E.fits')
except IOError:
    rn = '%.4f' % np.random.rand(1)[0]
    print('Fits files arleady exist.')
    print('New file name is '+rn)
    fits.PrimaryHDU(X, header).writeto(fits_name+'_X'+rn+'.fits')
    fits.PrimaryHDU(Y, header).writeto(fits_name+'_Y'+rn+'.fits')
    fits.PrimaryHDU(E, header).writeto(fits_name+'_E'+rn+'.fits') 
    
print(str(time.time()-start_timer)+'s')


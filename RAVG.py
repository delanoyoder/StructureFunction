#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 05:33:57 2019

@author: Delano Yoder
"""
# Timer
import time
start_timer = time.time()

import numpy as np
from astropy.io import fits

# File name
file_name = 'smc_atca+pks.imbin.fits'

# Importing fits header and data
fits_file = fits.open(file_name)
header, data = fits_file[0].header, fits_file[0].data

# Integrating and normalizing data into an image
while (len(data.shape) > 2):
    data = np.nansum(data, 0)
image = data/np.amax(data[np.where(~np.isnan(data))])

# Variables for image width and length
image_width, image_length = image.shape[0], image.shape[1]

# Kernel diameter and step size
kernel_diameter, step_size = 51, 1

# Number of possible kernel steps across the image and down the image
number_x_kernels = int((image_width-kernel_diameter)/step_size+1)
number_y_kernels = int((image_length-kernel_diameter)/step_size+1)

# Allocating memory for the results
A = np.empty((number_x_kernels, number_y_kernels))*np.nan

# Variable for time estimate
percentage_counter = 100

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
            
            # Averaging the kernel
            a = np.nanmean(kernel)

            # Storing the results
            A[x, y] = a
            
        # Time estimate
        if (x == 1 and y == 0):
            time_estimate = time.time()-start_timer
            time_scaling_factor = number_x_kernels*number_y_kernels/(number_y_kernels+1)
            time_estimate = time_estimate*time_scaling_factor
            if (time_estimate  <= 60):
                print('Estimated time is '+str(time_estimate)+'secs')
            elif (time_estimate > 60 and time_estimate <= 3600):
                print('Estimated time is '+str(time_estimate/60.)+'mins')
            else:
                print('Estimated time is '+str(time_estimate/3600.)+'hrs')
            percentage_counter = 0
        
        # Percentage complete display
        percentage_done = 100.*(x*number_y_kernels+y+1)/(number_x_kernels*number_y_kernels-1)
        if ( percentage_done >= percentage_counter):      
            print(str(int(100.*(x*number_y_kernels+y+1)/(number_x_kernels*number_y_kernels-1)))+'%')
            percentage_counter += 10
                
            
# Variable for new fits file name
try:
    fits_name = header['OBJECT'][:3].upper()+'_'+header['TELESCOP'].upper()
    fits_name += '_INT+RAVG_'+str(kernel_diameter)+'_'+str(step_size)
except KeyError:
    fits_name = file_name[:-5]+'_INT+RAVG_'+str(kernel_diameter)+'_'+str(step_size)
    
# Creating the new fits file containing the structure function results
try:
    fits.PrimaryHDU(A, header).writeto(fits_name+'.fits')
except IOError:
    rn = '%.4f' % np.random.rand(1)[0]
    print('Fits files arleady exist.')
    print('New file name is '+rn)
    fits.PrimaryHDU(A, header).writeto(fits_name+rn+'.fits')
    
print(str(time.time()-start_timer)+'s')
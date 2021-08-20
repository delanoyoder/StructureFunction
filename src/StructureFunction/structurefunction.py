"""
@author: Delano Yoder
"""

import time
import math
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as pltc

class Fits:

    def __init__(self, file_name):
        self.file_name = file_name
        self.get_header()
        self.get_data()

    def get_header(self):
        with fits.open(self.file_name) as ff:
            self.header = ff[0].header

    def get_data(self):
        with fits.open(self.file_name) as ff:
            self.data = ff[0].data

    def get_image(self):
        if self.data.ndim == 3:
            image = np.nansum(self.data, 0)
        else:
            image = self.data
        
        non_nan = np.where(~np.isnan(image))
        self.image = image / np.amax(image[non_nan])

    def get_structure_function(self, num_bins=1):
        self.structure_function = SF(self.data, num_bins)

    def get_rolling_structure_function(self):
        pass

    def get_SCA_structure_function(self):
        pass

    def get_SCA_rolling_structure_function(self):
        pass

class Image:

    def __init__(self, file_name):
        self.image = data

class SF:

    def __init__(self, data, num_bins):
        self.data = data
        self.get_image()
        self.get_attributes()
        self.get_sf()
        if num_bins > 1:
            self.bin(num_bins)

    def get_image(self):
        if self.data.ndim == 3:
            image = np.nansum(self.data, 0)
        else:
            image = self.data
        
        non_nan = np.where(~np.isnan(image))
        self.image = image / np.amax(image[non_nan])

    def get_attributes(self):
        self.width, self.length = self.image.shape
        self.max_distance = min(self.width / 2.0, self.length / 2.0)

    def get_sf(self):
        self.dict = {}
        self.parallels()
        self.diagonals()
        self.sort()

    def parallels(self):
        w = self.width
        l = self.length
        md = self.max_distance

        for d in range(1, int(md)):

            if d <= md:
                
                right = np.nanmean((self.image[0:w-d,0:l] - self.image[d:w,0:l])**2)
                up = np.nanmean((self.image[0:w,0:l-d] - self.image[0:w,d:l])**2)

                if (d in self.dict):
                    self.dict[d].append(right)
                    self.dict[d].append(up)
                else:
                    self.dict[d] = [right, up]
    
    def diagonals(self):
        w = self.width
        l = self.length
        md = self.max_distance

        for x in range(1, int(md)):
            for y in range(1, int(md)):
                
                d = (x**2+y**2)**0.5
                
                if d <= md:
                    
                    up_right = np.nanmean((self.image[0:w-x,0:l-y] - self.image[x:w,y:l])**2)
                    up_left = np.nanmean((self.image[x:w,0:l-y] - self.image[0:w-x,y:l])**2)
                
                    if (d in self.dict):
                        self.dict[d].append(up_right)
                        self.dict[d].append(up_left)
                    else:
                        self.dict[d] = [up_right, up_left]

    def sort(self):
        distances = []
        values = []
        errors = []
        
        for d in sorted(self.dict.keys()):
            distances.append(d)
            values.append(np.nanmean(self.dict[d]))
            errors.append(np.nanstd(self.dict[d]))

        self.distances = np.array(distances)
        self.values = np.array(values)
        self.errors = np.array(errors)
            


    def bin(self, num_bins):
        self.binned_distances = []
        self.binned_values = []
        self.binned_errors = []
        
        bin_spacing = np.logspace(0, max(self.distances)+1, num=num_bins+1)
        for i in range(num_bins):
            in_bin = (self.distances >= bin_spacing[i]) * (self.distances < bin_spacing[i])
            self.binned_distances.append(np.nanmean(self.distances[in_bin]))
            print(np.nansum(self.values[in_bin] / self.errors[in_bin]**2))
            self.binned_values.append(np.nansum(self.values[in_bin] / self.errors[in_bin]**2) 
                                    / np.nansum(1 / self.errors[in_bin]**2))
            self.binned_errors.append((1 / np.nansum(1 / self.errors[in_bin]**2))**0.5)

class RSF:

    def __init__(self):
        pass

class SCASF:
    def __init__(self):
        pass

class SCARSF:
    def __init__(self):
        pass
"""
@author: Delano Yoder
"""

import time
import math
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as pltc

# This Fits class opens and stores the contents of a .fits file. This class has 
# methods that call other objects to perform structure function analysis on.
# This includes SF, RSF, SCASF, and SCARSF. 
class Fits:

    def __init__(self, file_name):
        self.file_name = file_name
        self.get_header()
        self.get_data()

    # Stores FITS header.
    def get_header(self):
        with fits.open(self.file_name) as ff:
            self.header = ff[0].header

    # Stores FITS data.
    def get_data(self):
        with fits.open(self.file_name) as ff:
            self.data = ff[0].data

    # Stores the normalized integrated image of the FITS data.
    def get_image(self):

        # Integrate image if there are velocity channels.
        if self.data.ndim == 3:
            image = np.nansum(self.data, 0)
        else:
            image = self.data
        
        # Normalize the image.
        non_nan = np.where(~np.isnan(image))
        self.image = image / np.amax(image[non_nan])

    # Performs the structure function of the whole image. Keyword argument 
    # num_bins will bin the structure function distances logarithmically into
    # the specified amount of bins. Keyword argument max_distance will cut of
    # structure function distances larger than the given distance.
    def get_structure_function(self, num_bins=1, max_distance=None):
        self.get_image()
        self.structure_function = SF(self.image, num_bins, max_distance)

    # To do
    def get_rolling_structure_function(self):
        pass

    # To do
    def get_SCA_structure_function(self):
        pass

    # To do
    def get_SCA_rolling_structure_function(self):
        pass

# This Image class opens and stores the contents of a .fits file. This class has 
# methods that call other objects to perform structure function analysis on.
# This includes SF, RSF, SCASF, and SCARSF. 
class Image:

    def __init__(self, file_name):
        self.image = file_name
        self.normalize_image()

    # Normalizes the image.
    def normalize_image(self):
        non_nan = np.where(~np.isnan(self.image))
        self.image = self.image / np.amax(self.image[non_nan])

    # Performs the structure function of the whole image. Keyword argument 
    # num_bins will bin the structure function distances logarithmically into
    # the specified amount of bins. Keyword argument max_distance will cut of
    # structure function distances larger than the given distance.
    def get_structure_function(self, num_bins=1, max_distance=None):
        self.structure_function = SF(self.image, num_bins, max_distance)

# This class takes in either a Fits or Image object and uses its data to perform
# and return a structure function analysis on the whole image.
class SF:

    def __init__(self, image, num_bins, max_distance):
        self.image = image
        self.get_sf()

        # Bin the structure function distances, values, and errors if the user
        # passes through a valid argument.
        if num_bins > 1:
            self.bin(num_bins, max_distance)

    # Structure function procedure.
    def get_sf(self):
        self.dict = {}
        self.right()
        self.up()
        self.diagonals()
        self.sort()

    # Calculating the squared differences on horizontal pixel separations.
    def right(self):
        w, l = self.image.shape

        for x in range(1, w):
                
            right = np.nanmean((self.image[0:w-x,0:l] - self.image[x:w,0:l])**2)

            if (x in self.dict):
                self.dict[x].append(right)
            else:
                self.dict[x] = [right]

    # Calculating the squared differences on vertical pixel separations.
    def up(self):
        w, l = self.image.shape

        for y in range(1, l):
                
            up = np.nanmean((self.image[0:w,0:l-y] - self.image[0:w,y:l])**2)

            if (y in self.dict):
                self.dict[y].append(up)
            else:
                self.dict[y] = [up]
    
    # Calculating the squared differences on diagonal pixel separations.
    def diagonals(self):
        w, l = self.image.shape

        for x in range(1, w):
            for y in range(1, l):
                
                d = (x**2+y**2)**0.5
                    
                up_right = np.nanmean((self.image[0:w-x,0:l-y] - self.image[x:w,y:l])**2)
                up_left = np.nanmean((self.image[x:w,0:l-y] - self.image[0:w-x,y:l])**2)
            
                if (d in self.dict):
                    self.dict[d].append(up_right)
                    self.dict[d].append(up_left)
                else:
                    self.dict[d] = [up_right, up_left]

    # Sorting all structure function data by distance.
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

    def bin(self, num_bins, max_distance):
        self.binned_distances = []
        self.binned_values = []
        self.binned_errors = []
        
        if max_distance == None:
            max_distance = max(self.distances)

        bin_spacing = np.logspace(0, np.log10(max_distance), num=num_bins+1)
        bin_spacing[-1] += 1
        for i in range(num_bins):
            in_bin = (self.distances >= bin_spacing[i]) * (self.distances < bin_spacing[i+1])
            if True in in_bin:
                self.binned_distances.append(np.nanmean(self.distances[in_bin]))
                self.binned_values.append(np.nansum(self.values[in_bin] / self.errors[in_bin]**2) 
                                        / np.nansum(1 / self.errors[in_bin]**2))
                self.binned_errors.append((1 / np.nansum(1 / self.errors[in_bin]**2))**0.5)

# To do
class RSF:

    def __init__(self):
        pass

# To do
class SCASF:
    def __init__(self):
        pass

# To do
class SCARSF:
    def __init__(self):
        pass
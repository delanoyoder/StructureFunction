"""
@author: Delano Yoder
"""

import time
import math
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as pltc

class Image:

    def __init__(self, file_name):
        self.file_name = file_name

        # Reading fits file and getting fits header and data.
        with fits.open(file_name) as ff:
            self.header = ff[0].header
            self.data = ff[0].data

        # Place holders for structure function attributes.
        self.structure_function = None
        self.rolling_structure_function = None
        self.SCA_structure_function = None
        self.SCA_rolling_structure_function = None

    def get_structure_function(self, num_bins=0):
        self.integrate_image()
        self.image_attributes()
        self.parallels()
        self.diagonals()

        if num_bins == 0:
            self.structure_function = SF(self.SF_dict)
        else:
            self.structure_function = SF(self.SF_dict).bin(num_bins)

    def get_rolling_structure_function(self):
        pass

    def get_SCA_structure_function(self):
        pass

    def get_SCA_rolling_structure_function(self):
        pass

    def integrate_image(self):
        # Integrating if data is a PPV cube.
        if self.data.ndim == 3:
            image = np.nansum(self.data, 0)
        else:
            image = self.data
            
        # Normalizing image, ignoring NaN values.
        non_nan = np.where(~np.isnan(image))
        self.image = image / np.amax(image[non_nan])

    def image_attributes(self):
        self.SF_dict = {}

        self.width, self.length = self.image.shape
        self.max_distance = min(self.width / 2.0, self.length / 2.0)

    def parallels(self):
        w = self.width
        l = self.length
        md = self.max_distance

        for d in range(1, int(md)):

            if d <= md:
                
                right = np.nanmean((self.image[0:w-d,0:l] - self.image[d:w,0:l])**2)
                up = np.nanmean((self.image[0:w,0:l-d] - self.image[0:w,d:l])**2)

                if (d in self.SF_dict):
                    self.SF_dict[d].append(right)
                    self.SF_dict[d].append(up)
                else:
                    self.SF_dict[d] = [right, up]
    
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
                
                    if (d in self.SF_dict):
                        self.SF_dict[d].append(up_right)
                        self.SF_dict[d].append(up_left)
                    else:
                        self.SF_dict[d] = [up_right, up_left]

class SF:

    def __init__(self, SF_dict):
        self.distances = []
        self.values = []
        self.errors = []
        
        for d in sorted(SF_dict.keys()):
            self.distances.append(d)
            self.values.append(np.nanmean(SF_dict[d]))
            self.errors.append(np.nanstd(SF_dict[d]))

    def bin(self, num_bins):
        self.binned_distances = []
        self.binned_values = []
        self.binned_errors = []
        
        bin_spacing = np.logspace(0, max(self.distances)+1, num=num_bins+1)
        for i in range(num_bins):
            in_bin = (self.distances >= bin_spacing[i]) * (self.distances < bin_spacing[i])
            self.binned_distances.append(np.nanmean(self.distances[in_bin]))
            self.binned_values.append(np.nansum(self.values[in_bin] / self.errors[in_bin]**2) 
                                    / np.nansum(1 / self.errors[in_bin]**2))
            self.binned_errors.append((1 / np.nansum(1 / np.errors[in_bin]**2))**0.5)

def plot(image):

    if type(image) == Image:
        if ~hasattr(image, 'image'):
            image.integrate_image()
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot((111))
        if np.nanmean(image.image)/np.max(image.image) > 0.01:
            im = ax.imshow(image.image)
        else:
            im = ax.imshow(image.image, norm=pltc.LogNorm())
        if 'OBJECT' in image.header:
            ax.set_title(image.header['OBJECT'])
        if 'TARGNAME' in image.header:
            ax.set_title(image.header['TARGNAME'])
        ax.invert_yaxis()
        bbar = plt.colorbar(im)
        plt.show()


    if type(image) == SF:
        plt.figure(figsize=(8, 8))
        plt.axes(xscale='log')
        plt.xlabel('Linear Scale [pc]')
        plt.ylabel("Log(SF(r'))")
        if hasattr(image, 'binned_distances:'):
            plt.plot(image.binned_distances, np.log10(image.binned_values), 'bs', markersize=7, markerfacecolor='None')
        else:
            plt.plot(image.distances, np.log10(image.values), 'b.', markersize=5, markerfacecolor='None')
        plt.show()
            
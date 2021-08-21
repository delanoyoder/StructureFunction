"""
@author: Delano Yoder
"""

import time
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

        # Stepping through each horizontal pixel separation.
        for x in range(1, w):
                
            # Calculating the average of the squared difference of all pixel 
            # separations at this distance.
            right = np.nanmean((self.image[0:w-x,0:l] - self.image[x:w,0:l])**2)

            # Appending the average to the dictionary using the distance of 
            # the pixel separation as the key.
            if (x in self.dict):
                self.dict[x].append(right)
            else:
                self.dict[x] = [right]

    # Calculating the squared differences on vertical pixel separations.
    def up(self):
        w, l = self.image.shape

        # Stepping through each vertical pixel separation.
        for y in range(1, l):
                
            # Calculating the average of the squared difference of all pixel 
            # separations at this distance.
            up = np.nanmean((self.image[0:w,0:l-y] - self.image[0:w,y:l])**2)

            # Appending the average to the dictionary using the distance of 
            # the pixel separation as the key.
            if (y in self.dict):
                self.dict[y].append(up)
            else:
                self.dict[y] = [up]
    
    # Calculating the squared differences on diagonal pixel separations.
    def diagonals(self):
        w, l = self.image.shape

        # Stepping through each diagonal pixel separation.
        for x in range(1, w):
            for y in range(1, l):
                
                # Calculating the distance of the pixel separation.
                d = (x**2+y**2)**0.5
                
                # Calculating the average of the squared difference of all pixel 
                # separations at this distance.
                up_right = np.nanmean((self.image[0:w-x,0:l-y] - self.image[x:w,y:l])**2)
                up_left = np.nanmean((self.image[x:w,0:l-y] - self.image[0:w-x,y:l])**2)
            
                # Appending the average to the dictionary using the distance of 
                # the pixel separation as the key.
                if (d in self.dict):
                    self.dict[d].append(up_right)
                    self.dict[d].append(up_left)
                else:
                    self.dict[d] = [up_right, up_left]

    # Sorting all structure function data by distance.
    def sort(self):

        # Getting number of discrete distances.
        n = len(self.dict)

        # Allocating memory for structure function data arrays.
        self.distances = np.empty((0,n), float)
        self.values = np.empty((0,n), float)
        self.errors = np.empty((0,n), float)
        
        # Storing structure function data arrays sorted by distance.
        for d in sorted(self.dict.keys()):
            self.distances = np.append(self.distances, np.array(d))
            self.values = np.append(self.values, np.array(np.nanmean(self.dict[d])))
            self.errors = np.append(self.errors, np.array(np.nanstd(self.dict[d])))

    # Binning structure function data logarithmically based on number of bins, 
    # given by num_bins argument, and the maximum distance, given by max_distance.
    # A maximum distance might be needed since the structure function becomes 
    # non-linear at higher pixel separations.
    def bin(self, num_bins, max_distance):

        # Empty lists for binned structure function data.
        self.binned_distances = []
        self.binned_values = []
        self.binned_errors = []
        
        # If a max_distance argument isn't given use all distances.
        if max_distance == None:
            max_distance = max(self.distances)

        # Calculate the logarithmic spacing of the distances.
        bin_spacing = np.logspace(0, np.log10(max_distance), num=num_bins+1)
        bin_spacing[-1] += 1 # Allows for the last distance to be included in the last bin.

        # Step through each bin and store the structure function data in that bin.
        for i in range(num_bins):

            # Getting truth table for the indicies of the structure function data in the bin.
            in_bin = (self.distances >= bin_spacing[i]) * (self.distances < bin_spacing[i+1])

            # Checking that there is at least one data point in the bin and storing structure
            # function data if so.
            if True in in_bin:
                self.binned_distances.append(np.nanmean(self.distances[in_bin]))

                # Statistically weighting each value in the bin based on that distance's standard deviation.
                self.binned_values.append(np.nansum(self.values[in_bin] / self.errors[in_bin]**2) 
                                        / np.nansum(1 / self.errors[in_bin]**2))
                    
                # Statistically binning the standard deviation data.
                self.binned_errors.append((1 / np.nansum(1 / self.errors[in_bin]**2))**0.5)

# To do
class RSF:

    def __init__(self, image, kernel_radius, step_size, num_bins, max_distance):
        self.image = image
        self.kernel = kernel_radius
        self.step = step_size
        self.get_rsf()

    # Rolling structure function procedure.
    def get_sf(self):
        self.sample()
        self.roll_kernel()
        
        self.sort()

    def sample(self):
        pass

    def roll_kernel(self):

        w, l = self.image.shape

        # Number of kernels along the x-axis and y-axis of the image.
        num_x_kernels = int(((w - (self.kernel*2+1)) / self.step) + 1)
        num_y_kernels = int(((l - (self.kernel*2+1)) / self.step) + 1)

                # Stepping through each kernel
        for x in range(num_x_kernels):
            for y in range(num_y_kernels):
                
                # Variables for the position of each edge of the kernel
                top = y * self.step
                left = x * self.step
                bottom = (y * self.step) + (self.kernel*2+1)
                right = (x * self.step) + (self.kernel*2+1)
                
                # Checking if kernel is within the image
                if (right <= w and bottom <= l):
                    kernel = image[left:right, top:bottom].copy()
                    
                    # Calculating maximum limit
                    max_r = (self.kernel*2+1)/2.0
                    
                    # Making the kernal circlular   
                    r_2 = (np.arange(0, (self.kernel*2+1), dtype=float)-max_r)**2
                    cdist = np.empty(((self.kernel*2+1), (self.kernel*2+1)))
                    for i in range((self.kernel*2+1)):
                        cdist[0:, i] = (r_2+r_2[i])**0.5
                    kernel[cdist > max_r] = np.nan
                    
                    # Performing structure function on the kernel
                    x_SF, y_SF, e_SF = SF(kernel)
                    
                    # Allocating memory for the results
                    if (x+y == 0):
                        X = np.empty((len(y_SF)))*np.nan
                        Y = np.empty((num_x_kernels, num_y_kernels, len(y_SF)))*np.nan
                        E = np.empty((num_x_kernels, num_y_kernels, len(e_SF)))*np.nan
                    
                    # Storing the results
                    X = x_SF
                    Y[x, y, :] = y_SF
                    E[x, y, :] = e_SF

    def get_kernel(self, x, y, circular=False):
        pass



# To do
class SCASF:
    def __init__(self):
        pass

# To do
class SCARSF:
    def __init__(self):
        pass
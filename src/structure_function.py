"""
@author: Delano Yoder
"""

import numpy as np
from skimage.util.shape import view_as_windows


# This class takes in either a Fits or Image object and uses its data to perform
# and return a structure function analysis on the whole image.
class SF:
    def __init__(self, image, num_bins, max_distance=None):
        self.image = self.normalize_image(image)
        self.get_sf()

        # Bin the structure function distances, values, and errors if the user
        # passes through a valid argument.
        if num_bins is not None:
            self.bin(num_bins, max_distance)

        self.get_slope()

    @staticmethod
    def normalize_image(image):
        return image / np.nanmax(image)

    # Structure function procedure.
    def get_sf(self):
        self.dict = {}
        self.right()
        self.up()
        self.diagonals()
        self.sort()

    def get_slope(self):
        # Convert the data to log scale
        log_binned_distances = np.log10(self.binned_distances)
        log_binned_values = np.log10(self.binned_values)

        # Use polyfit with degree 1 to perform a linear regression
        self.slope, self.intercept = np.polyfit(
            log_binned_distances, log_binned_values, 1
        )

    # Calculating the squared differences on horizontal pixel separations.
    def right(self):
        w, l = self.image.shape

        # Stepping through each horizontal pixel separation.
        for x in range(1, w):
            # Calculating the average of the squared difference of all pixel
            # separations at this distance.
            right = np.nanmean((self.image[0 : w - x, 0:l] - self.image[x:w, 0:l]) ** 2)

            # Appending the average to the dictionary using the distance of
            # the pixel separation as the key.
            if x in self.dict:
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
            up = np.nanmean((self.image[0:w, 0 : l - y] - self.image[0:w, y:l]) ** 2)

            # Appending the average to the dictionary using the distance of
            # the pixel separation as the key.
            if y in self.dict:
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
                d = (x**2 + y**2) ** 0.5

                # Calculating the average of the squared difference of all pixel
                # separations at this distance.
                up_right = np.nanmean(
                    (self.image[0 : w - x, 0 : l - y] - self.image[x:w, y:l]) ** 2
                )
                up_left = np.nanmean(
                    (self.image[x:w, 0 : l - y] - self.image[0 : w - x, y:l]) ** 2
                )

                # Appending the average to the dictionary using the distance of
                # the pixel separation as the key.
                if d in self.dict:
                    self.dict[d].append(up_right)
                    self.dict[d].append(up_left)
                else:
                    self.dict[d] = [up_right, up_left]

    # Sorting all structure function data by distance.
    def sort(self):
        # Getting number of discrete distances.
        n = len(self.dict)

        # Allocating memory for structure function data arrays.
        self.distances = np.empty((0, n), float)
        self.values = np.empty((0, n), float)
        self.errors = np.empty((0, n), float)

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
            max_distance = np.power(10, np.log10(max(self.distances)) / 2)

        # Calculate the logarithmic spacing of the distances.
        bin_spacing = np.logspace(0, np.log10(max_distance), num=num_bins + 1)
        # Allows for the last distance to be included in the last bin.
        bin_spacing[-1] += 1

        # Step through each bin and store the structure function data in that bin.
        for i in range(num_bins):
            # Getting truth table for the indicies of the structure function data in the bin.
            in_bin = (self.distances >= bin_spacing[i]) * (
                self.distances < bin_spacing[i + 1]
            )

            # Checking that there is at least one data point in the bin and storing structure
            # function data if so.
            if True in in_bin:
                self.binned_distances.append(np.nanmean(self.distances[in_bin]))

                # Statistically weighting each value in the bin based on that distance's standard deviation.
                self.binned_values.append(
                    np.nansum(self.values[in_bin] / self.errors[in_bin] ** 2)
                    / np.nansum(1 / self.errors[in_bin] ** 2)
                )

                # Statistically binning the standard deviation data.
                self.binned_errors.append(
                    (1 / np.nansum(1 / self.errors[in_bin] ** 2)) ** 0.5
                )


# To do
class RSF:
    def __init__(self, image, kernel_size, step_size, num_bins, max_distance=None):
        self.get_rsf(image, kernel_size, step_size, num_bins)

    def get_rsf(self, image, k, s, num_bins):
        # Diameter will be used for window shape
        kernel_diameter = 2 * k + 1

        # Create overlapping windows on the image
        kernels = view_as_windows(image, (kernel_diameter, kernel_diameter), step=s)

        # Initialize result array
        self.results = np.zeros(kernels.shape[:2])

        # iterate over each window
        for i in range(kernels.shape[0]):
            for j in range(kernels.shape[1]):
                # Extract window
                kernel = kernels[i, j, :, :]

                # Check if window is a circular kernel
                y, x = np.ogrid[-k : k + 1, -k : k + 1]
                mask = x**2 + y**2 <= k**2

                # Convert mask's False values to NaN
                mask = mask.astype(float)
                mask[mask == 0] = np.nan

                # Apply circular mask to the window
                kernel = kernel * mask

                # Apply get_sf function to kernel and store result
                self.results[i, j] = SF(kernel, num_bins).slope
                print(f"({i}, {j}): {self.results[i, j]}")

        import matplotlib.pyplot as plt

        plt.imshow(image)
        plt.show()
        plt.imshow(self.results)
        plt.show()
        breakpoint()


# To do
class SCASF:
    def __init__(self):
        pass


# To do
class SCARSF:
    def __init__(self):
        pass

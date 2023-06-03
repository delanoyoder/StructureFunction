from astropy.io import fits
import numpy as np


class Fits:
    """This Fits class opens and stores the contents of a .fits file. This class
    has methods that call other objects to perform structure function analysis
    on. This includes SF, RSF, SCASF, and SCARSF.
    """

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

import matplotlib.pyplot as plt
from structurefunction import *

class Plot:

    def __init__(self, fits):
        self.fits = fits

    def plot_image(self):

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot((111))
        im = ax.imshow(self.fits.image)
        if 'OBJECT' in self.fits.header:
            ax.set_title(self.fits.header['OBJECT'])
        if 'TARGNAME' in self.fits.header:
            ax.set_title(self.fits.header['TARGNAME'])
        ax.invert_yaxis()
        bbar = plt.colorbar(im)
        plt.show()



    def plot_sf(self):
        plt.figure(figsize=(8, 8))
        plt.axes(xscale='log')
        plt.xlabel('Linear Scale [pixels]')
        plt.ylabel("Log(SF(r'))")
        if hasattr(self.fits.structure_function, 'binned_distances'):
            plt.plot(self.fits.structure_function.binned_distances, 
                    np.log10(self.fits.structure_function.binned_values), 
                    'bs', markersize=7, markerfacecolor='None')
        else:
            plt.plot(self.fits.structure_function.distances, 
                    np.log10(self.fits.structure_function.values), 
                    'b.', markersize=5, markerfacecolor='None')
        plt.show()
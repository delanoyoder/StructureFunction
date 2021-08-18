import matplotlib.pyplot as plt
from structurefunction import *

class Plot:

    def __init__(self, image):
        self.image = image
        
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
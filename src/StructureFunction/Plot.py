import matplotlib.pyplot as plt
from structurefunction import *

class Plot:

    def __init__(self, image):
        self.image = image

    def plot(self):
        
        if type(self.image) == Image:
            if ~hasattr(self.image, 'image'):
                self.image.integrate_image()
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot((111))
            if np.nanmean(self.image.image)/np.max(self.image.image) > 0.01:
                im = ax.imshow(self.image.image)
            else:
                im = ax.imshow(self.image.image, norm=pltc.LogNorm())
            if 'OBJECT' in self.image.header:
                ax.set_title(self.image.header['OBJECT'])
            if 'TARGNAME' in self.image.header:
                ax.set_title(self.image.header['TARGNAME'])
            ax.invert_yaxis()
            bbar = plt.colorbar(im)
            plt.show()


        if type(self.image) == SF:
            plt.figure(figsize=(8, 8))
            plt.axes(xscale='log')
            plt.xlabel('Linear Scale [pc]')
            plt.ylabel("Log(SF(r'))")
            if hasattr(self.image, 'binned_distances:'):
                plt.plot(self.image.binned_distances, np.log10(self.image.binned_values), 'bs', markersize=7, markerfacecolor='None')
            else:
                plt.plot(self.image.distances, np.log10(self.image.values), 'b.', markersize=5, markerfacecolor='None')
            plt.show()
import matplotlib.pyplot as plt
import numpy as np


class Plot:
    @staticmethod
    def plot_image(fits):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot((111))
        im = ax.imshow(fits.data)
        if "OBJECT" in fits.header:
            ax.set_title(fits.header["OBJECT"])
        if "TARGNAME" in fits.header:
            ax.set_title(fits.header["TARGNAME"])
        ax.invert_yaxis()
        bbar = plt.colorbar(im)
        plt.show()

    @staticmethod
    def plot_sf(sf):
        plt.figure(figsize=(8, 8))
        plt.axes(xscale="log")
        plt.xlabel("Linear Scale [pixels]")
        plt.ylabel("Log(SF(r'))")
        if hasattr(sf, "binned_distances"):
            plt.plot(
                sf.binned_distances,
                np.log10(sf.binned_values),
                "bs",
                markersize=7,
                markerfacecolor="None",
            )
        else:
            plt.plot(
                sf.distances,
                np.log10(sf.values),
                "b.",
                markersize=5,
                markerfacecolor="None",
            )
        plt.show()

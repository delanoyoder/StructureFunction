from structurefunction import *
from plot import *
import os

def test1(fits):
    fits.get_image()
    Plot(fits).plot_image()

def test2(fits):
    fits.get_structure_function()
    Plot(fits).plot_sf()

def test3(fits):
    fits.get_structure_function(50)
    Plot(fits).plot_sf()

def test4(fits):
    fits.get_structure_function(30, 6)
    Plot(fits).plot_sf()

root = "/Users/delanoyoder/Projects/StructureFunction/src/StructureFunction/data"
path = os.path.join(root, "sample.fits")
fits = Fits(path)
test1(fits)
test2(fits)
test3(fits)
test4(fits)
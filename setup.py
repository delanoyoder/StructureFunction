import os
import sys
import glob
import warnings
from setuptools import setup, Extension


file_required = "requirements.txt"

with open(file_required) as file_:
    required = file_.read().splitlines()

setup(
    name='StructureFunction',
    url='https://github.com/delanoyoder/StructureFunction.git',
    author='Delano Yoder',
    author_email='dayoder4@gmail.com',
    description='package for plotting a different structure function analysis on images',
    packages=['StructureFunction'],
    install_requires=required
)
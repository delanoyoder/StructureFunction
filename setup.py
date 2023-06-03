from setuptools import setup, find_packages

VERSION = "0.0.1"
DESCRIPTION = "Package for structure function applications"

# Read requirements.txt and store the packages into a list
with open("requirements.txt") as f:
    required_packages = f.read().splitlines()

# Setting up
setup(
    name="StructureFunction",
    version=VERSION,
    author="Delano Yoder",
    author_email="<dayoder4@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=required_packages,
    keywords=["python", "structure", "function", "fits", "analysis"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Image Processors",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
)

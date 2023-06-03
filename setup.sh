#!/bin/bash

# Create the virtual environment
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install your package using setup.py
pip install .

# Deactivate the virtual environment
deactivate

# If you want the virtual environment to remain active after the script, comment out the deactivate line.

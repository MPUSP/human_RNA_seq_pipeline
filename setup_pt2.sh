#!/bin/bash

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null
then
    echo "Python 3 could not be found. Please install Python 3 and try again."
    exit 1
fi

# Create a virtual environment in the ./venv directory
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Upgrade pip to the latest version
pip install --upgrade pip

# Install the required packages
pip install matplotlib scipy numpy

# Deactivate the virtual environment
deactivate

echo "Virtual environment setup complete. 'matplotlib', 'scipy', and numpy have been installed."

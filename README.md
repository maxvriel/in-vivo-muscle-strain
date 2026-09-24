[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21506263.svg)](https://doi.org/10.5281/zenodo.21506263)

# Time-resolved 3D imaging and strain analysis of the skeletal muscle

This repository contains the MATLAB code used to generate the figures in the paper:

**"Time-resolved 3D imaging and strain analysis of the skeletal muscle"**  
by Max H.C. van Riel, David G.J. Heesterbeek, Martijn Froeling, Tristan van Leeuwen, Cornelis A.T. van den Berg, and Alessandro Sbrizzi  
University Medical Center Utrecht, The Netherlands

## Requirements

- MATLAB with the Image Processing Toolbox (tested with Matlab 2023a)
- [ISMRMRD](https://github.com/ismrmrd/ismrmrd/)
- [elastix](https://elastix.dev/)

Make sure the paths in `setup.m` are configured correctly.

Note that there is a bug in the hdf5 library of Matlab 2025a and newer that can cause memory issues when loading multiple ISMRMRD files.

## Data

Download the dataset from [Zenodo](https://doi.org/10.5281/zenodo.17312307) and place it in a folder named `data`.

## Installation

1. Clone this repository.
2. Install the dependencies listed above.
3. Update the paths in `setup.m`.
4. Place the downloaded dataset in `data/`.

## Usage

To run the joint reconstruction, run the `runRecon.m` script.

Reconstruction may take several hours to complete.

After reconstruction, Figures 3–9 in the manuscript can be reproduced using the corresponding scripts.

## Registration outputs

For Figures 5–9, elastix registration outputs are stored in `registration/tmp`.

If the reconstructions or segmentation masks change, delete this folder before rerunning the figure generation scripts.

## Acknowledgements

This project uses the scientific colour maps **lajolla** and **vik**:

Crameri, F. (2018a). *Scientific colour maps*. Zenodo. https://doi.org/10.5281/zenodo.1243862

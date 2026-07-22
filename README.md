# Time-resolved 3D imaging and strain analysis for in vivo muscle dynamics
This repository includes the Matlab code used to produce the figures included in the paper titled "Time-resolved 3D imaging and strain analysis for in vivo muscle dynamics" by Max H.C. van Riel, David G.J. Heesterbeek, Martijn Froeling, Tristan van Leeuwen, Cornelis A.T. van den Berg, and Alessandro Sbrizzi from the University Medical Center Utrecht, The Netherlands.

## Installation
The Matlab code has been tested using Matlab 2023a with the Image Processing Toolbox.
To run the scripts, install the ISMRMRD code (https://github.com/ismrmrd/ismrmrd/) and elastix (https://elastix.dev/).
Make sure the paths in setup.m are set correctly.
Finally, download the dataset from https://doi.org/10.5281/zenodo.17312307 and place the dataset in a folder called 'data'.

## Usage
To run the joint reconstruction, run the script runRecon.m.
This will take several hours to complete.
After all datasets have been reconstructed, the figures 3-9 of the manuscript can be reproduced using the corresponding Matlab scripts.

## Acknowledgements
We make use of the scientific color maps lajolla and vik:
Crameri, F. (2018a), Scientific colour maps. Zenodo. doi.org/10.5281/zenodo.1243862.

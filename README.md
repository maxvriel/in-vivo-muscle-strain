[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21506263.svg)](https://doi.org/10.5281/zenodo.21506263)

# Time-resolved 3D imaging and strain analysis of the skeletal muscle

This repository contains the MATLAB code used to generate the figures in the paper:

**"Time-resolved 3D imaging and strain analysis of the skeletal muscle"**  
by Max H.C. van Riel, David G.J. Heesterbeek, Martijn Froeling, Tristan van Leeuwen, Cornelis A.T. van den Berg, and Alessandro Sbrizzi  
University Medical Center Utrecht, The Netherlands

## Requirements

- MATLAB 2023a with the Image Processing Toolbox
- [ISMRMRD](https://ismrmrd.readthedocs.io/en/latest/index.html) v1.14.3
- [elastix](https://elastix.dev/) 5.2.0

The code was tested with the indicated versions, but other versions may also work.
Note that there is a bug (confirmed by MathWorks support) in the hdf5 library of Matlab versions 2025a and newer (encountered using Matlab 2025b on Linux).
As a result, the script for Figure 5 may cause out-of-memory issues when run with affected Matlab versions because it loads multiple ISMRMRD files.

## Data

Download the dataset from [Zenodo](https://doi.org/10.5281/zenodo.17312307) and place it in a folder named `data`.

For the `runRecon.m` script and the scripts generating Figures 3, 4, 6, and 7, only the data from a single volunteer is required.
The scripts reproducing Figures 5 and 8 require all data to be downloaded and reconstructed.
The file `code.zip` from the dataset is not required for this repository, and the installation instructions at the end of the Zenodo page do not apply to this repository.

## Installation

1. Clone this repository and the ISMRMRD code:
    ```
    git clone https://github.com/ismrmrd/ismrmrd.git
    cd ismrmrd
    git checkout v1.14.3
    cd ..
    git clone https://github.com/maxvriel/in-vivo-muscle-strain.git
    ```
2. Add the `ismrmrd/matlab` folder to the Matlab path in `setup.m`.
3. Download the elastix release for your operating system from the [elastix releases page](https://github.com/SuperElastix/elastix/releases/tag/5.2.0).
4. Unzip the release and set `ELASTIXPATH` in `setup.m` to the `elastix-5.2.0` directory.
5. Download (part of) the dataset, unzip the folders and place them in `data/`.
6. Start MATLAB and run `runRecon.m` or (after reconstruction) one of the figure scripts.

## Usage

Check the recon settings at the top of `runRecon.m` and then run the script.
It solves the following optimization problem:

```math
\min_{\mathbf{m},\mathbf{v}_b}\mathcal{D}(\mathbf{m})+\lambda\mathcal{M}(\mathbf{m},\mathbf{v}_b)+\mu\mathcal{R}(\mathbf{v}_b)
```

with $\mathbf{m}$ the time-resolved image, $\mathbf{v}_b$ the spline coefficients of the velocity field, $\mathcal{D}$ the data consistency term, $\mathcal{M}$ the motion model, and $\mathcal{R}$ the regularization of the velocity field's spline coefficients.
It uses a multi-scale strategy, where first a reconstruction was performed at half the spatial resolution, which served as initialization for the reconstruction at the full spatial resolution.
More details about the reconstruction can be found in the paper.

Reconstruction may take several hours to complete, after which the results are saved to the `recon` folder.

After the required reconstructions have been performed, the figure scripts can be run. 
They reproduce the figures as shown in the paper.
- Figure 3: Reconstructed image slices and lines over time.
- Figure 4: Reconstructed images and displacement fields.
- Figure 5: Validation of the dynamic images and motion fields.
- Figure 6: Segmentation and octahedral shear strain maps.
- Figure 7: Mean octahedral shear strain per muscle, and muscle-wise differences between passive and contracted muscles.
- Figure 8: Mean octahedral shear strain differences per muscle for all scans.

## Registration outputs

For Figures 5–9, elastix registration outputs are stored in `registration/tmp`.

If the reconstructions or segmentation masks change, delete this folder before rerunning the figure generation scripts.

## Acknowledgements

This research was supported by the Dutch Research Council (NWO), grant 18897.

This project uses the scientific colormaps **lajolla** and **vik**:

Crameri, F. (2018a). *Scientific colour maps*. Zenodo. https://doi.org/10.5281/zenodo.1243862

The authors thank the ISMRM Reproducible Research Study Group for conducting a code review of the code (Version 1.0) supplied in this repository. The scope of the code review covered only the code's ease of download, quality of documentation, and ability to run, but did not consider scientific accuracy or code efficiency.

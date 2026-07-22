function writeRecon(fileName, images, timeStep, fieldOfView, vSplCoeffs, Bx, By, Bz)
%writeRecon Save reconstructed images and velocity fields to HDF5 file.
%
%   Inputs:
%   fileName: Path to the HDF5 file to save the reconstruction data
%   images: Reconstructed complex-valued images as a 4D array [x, y, z, t]
%   timeStep: Time step between image frames in seconds
%   fieldOfView: Field of view in mm
%   vSplCoeffs: Velocity spline coefficients as a 5D array 
%       [bx, by, bz, v, t]
%   Bx, By, Bz: Spline basis functions for each dimension, where each
%       column contains one B-spline
%   
% Copyright (c) 2026, UMC Utrecht
% Max van Riel, m.h.c.vanriel-3@umcutrecht.nl

complexTypeID = H5T.create('H5T_COMPOUND', 2*H5T.get_size('H5T_NATIVE_FLOAT'));
H5T.insert(complexTypeID, 'real', 0, 'H5T_NATIVE_FLOAT');
H5T.insert(complexTypeID, 'imag', H5T.get_size('H5T_NATIVE_FLOAT'), 'H5T_NATIVE_FLOAT');

fileID = H5F.create(fileName, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

reconID = H5G.create(fileID, 'recon', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');

[nx, ny, nz, nTime] = size(images, 1:4);
dataSpaceID = H5S.create_simple(4, [nTime, nz, ny, nx], []);
dataSetID = H5D.create(reconID, 'images', complexTypeID, dataSpaceID, 'H5P_DEFAULT');

images = struct('real', real(images), 'imag', imag(images));
H5D.write(dataSetID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', images);

h5create(fileName, '/recon/velocity/splines_x', size(Bx), 'Datatype', class(Bx))
h5write(fileName, '/recon/velocity/splines_x', Bx)
h5create(fileName, '/recon/velocity/splines_y', size(By), 'Datatype', class(By))
h5write(fileName, '/recon/velocity/splines_y', By)
h5create(fileName, '/recon/velocity/splines_z', size(Bz), 'Datatype', class(Bz))
h5write(fileName, '/recon/velocity/splines_z', Bz)
vSplCoeffs = reshape(vSplCoeffs, size(Bx, 2), size(By, 2), size(Bz, 2), [], nTime-1);
h5create(fileName, '/recon/velocity/coefficients', size(vSplCoeffs), 'Datatype', class(vSplCoeffs))
h5write(fileName, '/recon/velocity/coefficients', vSplCoeffs)

h5writeatt(fileName, '/recon', 'field_of_view', fieldOfView)
h5writeatt(fileName, '/recon', 'time_step', timeStep)

H5T.close(complexTypeID)
H5D.close(dataSetID)
H5S.close(dataSpaceID)
H5G.close(reconID)
H5F.close(fileID)

end

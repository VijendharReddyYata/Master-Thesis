function [ alpha, beta, gamma, R_yz, R_xz, R_xy, R ] = vol_orient_ACF( data, cube_size)
%VOL_ORIENT_CONVN extracts the orientation and scaling parameter of
% structures from given volume
% [ alpha beta gamma R] = vol_orient_convn( data ) extracts the orientation in
% three planes (yz, xz and xy) and scaling parameter
%
% 'alpha'   is a 3d matrix with orientation in YZ plane
% 'beta'    is a 3d matrix with orientation in xz plane
% 'gamma'   is a 3d matrix with orientation in XY plane 
% 'R_yz'       is a 3d matrix with scaling parameter in YZ plane
% 'R_xz'       is a 3d matrix with scaling parameter in xz plane
% 'R_xy'       is a 3d matrix with scaling parameter in XY plane
% 'R'       is a 3d matrix with scaling parameter (more anisotropic value / parameter among the 3 maximum projections of convolution matrix)
% Output arguments are 3d matrices with the data where exactly it is
% extracted from, which is reshaped according to the cube size given
%
% Input arguments
% 'data' should be a 3d volume
% 'cube_size' should be a parameter of positive value, which can be given
% according to the voxel size of the data
%
% It finds the number of rows, columns and slices of the input data and
% creates some preallocated matrices accordingly to store the results of
% alpha, beta, gamma and for anisotropy. Then, three for loops (slices, columns and rows) will be run throughout the
% volume with the given cube size. Function 'orient_convn' will be applied
% to every cube, so that orientation and anisotropy will be
% calculated and stores in the already preallocated matrices as a single
% array. Then the single arrays will be reshaped into a 3d volume according
% to the input data, from where the results came from. This will be done
% using the given 'cube size' in input.
%
%   Author:     Vijendhar Reddy Yata (vijaychinna.8352@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   17/05/2018
%   Last update:  26/05/2018

%% 1. For loop to extract fabric of cortical bone
A=size(data,1); % finding number of rows
B=size(data,2); % finding number of columns
C=size(data,3); % finding number of slices

yz_angle=zeros(1,ceil(A/cube_size*B/cube_size*C/cube_size)); % preallocation for alpha
xz_angle=zeros(1,ceil(A/cube_size*B/cube_size*C/cube_size)); % preallocation for beta
xy_angle=zeros(1,ceil(A/cube_size*B/cube_size*C/cube_size)); % preallocation for gamma
anisotropy_yz=zeros(1,ceil(A/cube_size*B/cube_size*C/cube_size)); % preallocation for R
anisotropy_xz=zeros(1,ceil(A/cube_size*B/cube_size*C/cube_size)); % preallocation for R
anisotropy_xy=zeros(1,ceil(A/cube_size*B/cube_size*C/cube_size)); % preallocation for R
anisotropy=zeros(1,ceil(A/cube_size*B/cube_size*C/cube_size)); % preallocation for R

count=1; % counter

for k=1:cube_size:C-cube_size % index for slices
    for j=1:cube_size:B-cube_size % index for columns
        for i=1:cube_size:A-cube_size % index for rows
            cube=data(i:i+cube_size-1, j:j+cube_size-1, k:k+cube_size-1); % loading cubes
            [alpha, beta, gamma, R_yz, R_xz, R_xy, R] = orient_ACF(cube); % to extract orientation and scaling property
            yz_angle(count)=alpha; % orientation in yz plane
            xz_angle(count)=beta; % orientation in xz plane
            xy_angle(count)=gamma; % orientation in xy plane
            anisotropy_yz(count)=R_yz;
            anisotropy_xz(count)=R_xz;
            anisotropy_xy(count)=R_xy;
            anisotropy(count)=R; % scaling property
            count=count+1; 
        end
    end
end

%% 2. Resizing and Reshaping of output matrices
%% 3.1. Resize
yz_angle=yz_angle(1,1:count-1); % resizing from a preallocated array 
xz_angle=xz_angle(1,1:count-1); % resizing from a preallocated array
xy_angle=xy_angle(1,1:count-1); % resizing from a preallocated array
anisotropy_yz=anisotropy_yz(1,1:count-1); % resizing from a preallocated array
anisotropy_xz=anisotropy_xz(1,1:count-1); % resizing from a preallocated array
anisotropy_xy=anisotropy_xy(1,1:count-1); % resizing from a preallocated array
anisotropy=anisotropy(1,1:count-1); % resizing from a preallocated array

%% 3.2. Reshape
rows=(A/cube_size); % scaling of data, depends on the 'cube_size'
rows_floor=floor(A/cube_size);
if rows==rows_floor
    rows=rows-1;
else
    rows=rows_floor;
end

columns=(B/cube_size); % scaling of data, depends on the 'cube_size'
columns_floor=floor(B/cube_size);
if columns==columns_floor
    columns=columns-1;
else
    columns=columns_floor;
end

slices=(C/cube_size); % scaling of data depends on the 'cube_size'
slices_floor=floor(C/cube_size); 
if slices==slices_floor
    slices=slices-1;
else
    slices=slices_floor;
end

alpha = reshape(yz_angle,[rows, columns, slices]); % Reshaping an array to a 3d volume
beta = reshape(xz_angle,[rows, columns, slices]); % Reshaping an array to a 3d volume
gamma = reshape(xy_angle,[rows, columns, slices]); % Reshaping an array to a 3d volume
R_yz = reshape(anisotropy_yz,[rows, columns, slices]); % Reshaping an array to a 3d volume
R_xz = reshape(anisotropy_xz,[rows, columns, slices]); % Reshaping an array to a 3d volume
R_xy = reshape(anisotropy_xy,[rows, columns, slices]); % Reshaping an array to a 3d volume
R = reshape(anisotropy,[rows, columns, slices]); % Reshaping an array to a 3d volume

end


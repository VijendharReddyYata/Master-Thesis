function [ alpha, beta, gamma, R_yz, R_xz, R_xy, R] = orient_ACF( data )
%ORIENT_CONVN extracts the orientation and scaling parameters of structures (orientation of low intensity pixels, which means pores in bone) from given data
%   [ alpha beta gamma R_yz R_xz R_xy R] = orient_convn( data ) extracts the orientation
%   from 3 planes(yz, xz and xy) of a 3D data and provides the scaling parameter (isotropic/anisotropic)
%
% 'alpha'   orientation in YZ plane
% 'beta'   orientation in xz plane
% 'gamma'   orientation in XY plane 
% 'R_yz'     anisotropy in yz plane
% 'R_xz'     anisotropy in xz plane
% 'R_xy'     anisotropy in xy plane
% 'R'       scaling parameter (more anisotropic value / parameter among the 3 maximum projections of convolution matrix)
% 
% Input data should be 3 dimentional and it can be of any class and size
%
% First, it inverts the intensities (low intensity voxels of pores becomes high intensity and voxels of bone becomes low intensity)
% If the input data is not in the class 'double', then data can be converted to 'double' class.
% To make sure all the elements orient around zero, 'sample scaling' will be
% done. Then smoothes the data using low pass filter. 
% Calcultes the convolution matrix and then maximum projections of convolution matrix in
% three planes (yz, xz and xy). Segmentation and binerization of the maximum intensity
% regions of maximum projections will be done. Using this data, orientation
% of region with  maximum area can be extracted from three different
% planes. Calculates the anisotropy (R) in three planes and shows the more
% anisotropic value in output along with the three planes anisotropy. 
% ______________________________________________________
% 
%   Author:     Vijendhar Reddy Yata (vijaychinna.8352@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   02/05/2018
%   Last update:  12/01/2018
% ______________________________________________________

%% 1. Image complement
data = double(data);
data = imcomplement(data);

%% 2. Denoising the data
%% 2.1. Smoothing data using low pass filter
LP=imgaussfilt3(data,1.2,'FilterSize',5); %sigma of 1.2 and kernel size 5*5*5 % Point spread function (PSF) 

%% 3. Sample scaling
% A = double(data); % converting class to double
[~, ~, ~]=size(data);
B=mean2(LP);
D=LP./B;
B1=mean2(D);
S=D-B1;
%% 2.2. Smoothing edges
%% 2.2.1. Creating a Gaussian kernel
% H=fspecial3('gaussian',[r c s]); % creates 3d gaussian kernel with the same size of input data
% 
% %% 2.2.1.1 Creating a Hamming kernel
% % H=hamming_3d(r, c, s); % creates 3d hamming kernel with the same size of input data
% 
% %% 2.2.2. Matrix multiplication of denoised data and Gaussian kernel
% MG=LP.*H;

%% 3. Calculate convolution matrix after matrix multiplication with GF
CC=convnfft(S,S(end:-1:1,end:-1:1,end:-1:1),'same'); % faster
% CC=convn(LP,LP(end:-1:1,end:-1:1,end:-1:1),'same'); % standard
max_CC = max(CC(:));
CCN = CC./max_CC; % normalize the data to 0 and 1
% CC=convnfft(MG,MG(end:-1:1,end:-1:1,end:-1:1),'same'); % faster

% tic
% CC=convn(MG,MG(end:-1:1,end:-1:1,end:-1:1),'same'); % slow
% toc

%% 4. Maximum projections of the convolution matrix
M1=permute((max(CCN,[],2)), [3 2 1]); %maximum projection of yz plane
M2=permute((max(CCN,[],1)), [1 3 2]); %maximum projection of xz plane
M3=max(CCN,[],3); %maximum projection of xy plane

%% 5. Segmentation of maximumm intensities
%% 5.1. Finding range for ROI
high=max(CCN(:)); %maximum value of ROI range

%% 5.2. Finding anlgle in yz plane (alpha)
alpha=NaN;
if ~isnan(M1) % if yz maximum projection is not NaN
%     thresh_yz = multithresh(squeeze(M1),1); % segmentation using Otsu's method
%     yz_low = min(thresh_yz); % minimum value of ROI range
%     [N,edges] = histcounts(M1);
%     CS = cumsum(N);
%     Y = prctile(CS,70);
%     [~, index] = min(abs(CS-Y));
%     yz_low = edges(1,index);
%     BW_yz=roicolor(squeeze(M1),yz_low, high);  
    BW_yz = imbinarize(squeeze(M1)); % segmentation of region with maximum intensity
    Rprop_yz = cell2mat(table2cell(regionprops('table',BW_yz,'Area','Orientation','MajorAxisLength','MinorAxisLength'))); % Characterizing regional properties
    [~,i] = max(Rprop_yz(:,1));
    alpha = Rprop_yz(i,4); % orientation angle in yz plane
end

%% 5.2. Finding anlgle in xz plane (beta)
beta=NaN;
if ~isnan(M2) % if xz maximum projection is not NaN
%     thresh_xz = multithresh(squeeze(M2),1);
%     xz_low = min(thresh_xz); % minimum value of ROI range
%     [N,edges] = histcounts(M2);
%     CS = cumsum(N);
%     Y = prctile(CS,70);
%     [~, index] = min(abs(CS-Y));
%     xz_low = edges(1,index);
%    BW_xz=roicolor(squeeze(M2),xz_low, high);  
    BW_xz = imbinarize(squeeze(M2)); % segmentation of region with maximum intensity
    Rprop_xz = cell2mat(table2cell(regionprops('table',BW_xz,'Area','Orientation','MajorAxisLength','MinorAxisLength'))); % Characterizing regional properties
    [~,j] = max(Rprop_xz(:,1));
    beta = Rprop_xz(j,4); % orientation angle in xz plane
end

%% 5.2. Finding anlgle in yz plane (alpha)
gamma=NaN;
if ~isnan(M3) % if xy maximum projection is not NaN
%     thresh_xy = multithresh(squeeze(M3),1);
%     xy_low = min(thresh_xy); % minimum value of ROI range
%     [N,edges] = histcounts(M3);
%     CS = cumsum(N);
%     Y = prctile(CS,70);
%     [~, index] = min(abs(CS-Y));
%     xy_low = edges(1,index);
%    BW_xy=roicolor(squeeze(M3),xy_low, high); 
    BW_xy = imbinarize(squeeze(M3)); % segmentation of region with maximum intensity 
    Rprop_xy = cell2mat(table2cell(regionprops('table',BW_xy,'Area','Orientation','MajorAxisLength','MinorAxisLength'))); % Characterizing regional properties
    [~,k] = max(Rprop_xy(:,1));
    gamma = Rprop_xy(k,4); % orientation angle in xy plane
end

angles = [alpha beta gamma]; % orientation in 3 planes
R_yz=NaN;
R_xz=NaN;
R_xy=NaN;
R=NaN;
if ~isnan(angles)
    R_yz=Rprop_yz(i,3)/Rprop_yz(i,2); % isotropic parameter in yz plane
    R_xz=Rprop_xz(j,3)/Rprop_xz(j,2); % isotropic parameter in xz plane
    R_xy=Rprop_xy(k,3)/Rprop_xy(k,2); % isotropic parameter in xy plane
    R_array=[R_yz R_xz R_xy]; % anisotropy of three planes
    R=min(R_array); % minimum isotropic parameter value (more anisotropic) among three planes
end
end

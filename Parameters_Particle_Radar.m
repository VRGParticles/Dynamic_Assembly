%% parameters for particle classification and tracking 

%% parameters for image_preprocess.m
%parameters for image preparation
%cropping
param.cropping = 1; %turn on and off cropping
param.crop_x1 = 0; %top left crop
param.crop_y1 = 0; %bottom left crop
param.crop_x2 = 1000; %top right crop
param.crop_y2 = 1000; %bottom right crop

%filter local standard deviation
%The particles are first segmented from the background by thresholding the
%particles based on the local standard deviation of the particles.
param.stdfilt_NHOOD = 11; %neighborhood for stdfilt function
param.stdfilt_threshold = 7; %Threshold for standard deviation filter

%Gaussian Filter
param.gaussian_sigma = 1; %standard deviation of Gaussian filter.

%Intensity Filter
param.image_adjust_LOW_IN =100; %Low range of the particles' intensities in the Gaussian-filtered image.
param.image_adjust_HIGH_IN =170; %High range of the particles' intensities in the Gaussian-filtered image.
param.image_adjust_LOW_OUT =0; %Low range of the particles' intensities in the newly scaled image (mask).
param.image_adjust_HIGH_OUT =255; %High range of the particles' intensities in the newly scaled image (mask).
param.image_intensity_threshold = 100; %threshold for the image intensity.

%90 180 0 255 105
%Particle Size Filter
param.area_filter_low = 70; %Smallest particle area
param.area_filter_high = 155; %Largest particle area
param.eccentricity_threshold = 0.7; %eccentricity threshold.
param.perimeter_threshold = 50; %perimeter threshold.

%% Parameters for angle differences
param.particle_distance_threshold = 17; %Threshold to determine the radius about each particle center where neighbors' centers fall.
param.malmin = 0; %minimum length a particle must be (in major axis), set to 0 if you don't want this
param.areamin = 0; %minimum a_c_area a particle must have, set to 0 if you don't want this
param.perimin = 0; %minimum perimeter a particle must have. set to 0 if you don't want this
param.e2c =  100; %cutoff distance of ends of each particle to the center of a neighbor.
param.max_chain_length = 25; %the maximum chain length expected

%% Parameters for Parallel Angles
%for the chains, we consider parallel angles.
param.para_angle1_low = 0; %these are given as examples, will not be true if absolute value of the angles are considered.
param.para_angle1_high = 30;

%second set of angles for set 2.
param.para_angle2_low = 150;
param.para_angle2_high = 180;

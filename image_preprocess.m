function a_i_image_eccen = image_preprocess(a_i_Im_original, param)

% ; NAME:
% ; image_preprocess
% ; PURPOSE:
% ; Convert images or consecutive video frames to binary images - particles
% ; background.
% ; CATEGORY:
% ; Image Processing
% ; CALLING SEQUENCE:
% ; a_i_image_eccen = image_preprocess(a_i_Im_original, param)
% ; INPUTS:
% ; a_i_Im_original: the unedited image / frame of the video
%   param:  a structure containing parameters for the various parts of the
%   image editing process: cropping, filtering, etc. Each is defined in the comments of the parameter
%   file.
% ; OUTPUTS:
% ; imageout: the binary image.


 %rgb check
 
 %Parameter read in

 %parameters for image preparation
%cropping
cropping = param.cropping;  %turn on and off cropping
crop_x1 = param.crop_x1; %top left crop
crop_y1 = param.crop_y1; %bottom left crop
crop_x2 = param.crop_x2; %top right crop
crop_y2 = param.crop_y2; %bottom right crop

%filter local standard deviation
%The particles are first segmented from the background by thresholding the
%particles based on the local standard deviation of the particles.
stdfilt_NHOOD = param.stdfilt_NHOOD; %neighborhood for stdfilt function
stdfilt_threshold = param.stdfilt_threshold; %Threshold for standard deviation filter

%Gaussian Filter
gaussian_sigma = param.gaussian_sigma; %standard deviation of Gaussian filter.

%Intensity Filter
image_adjust_LOW_IN = param.image_adjust_LOW_IN; %Low range of the particles' intensities in the Gaussian-filtered image.
image_adjust_HIGH_IN = param.image_adjust_HIGH_IN; %High range of the particles' intensities in the Gaussian-filtered image.
image_adjust_LOW_OUT = param.image_adjust_LOW_OUT; %Low range of the particles' intensities in the newly scaled image (mask).
image_adjust_HIGH_OUT = param.image_adjust_HIGH_OUT; %High range of the particles' intensities in the newly scaled image (mask).
image_intensity_threshold = param.image_intensity_threshold; %threshold for the image intensity.

%Particle Size Filter
area_filter_low = param.area_filter_low; %Smallest particle area
area_filter_high = param.area_filter_high; %Largest particle area
eccentricity_threshold = param.eccentricity_threshold ; %eccentricity threshold.
perimeter_threshold = param.perimeter_threshold ; %perimeter threshold.
 
% Cropping and color check
    rgb_check = size(a_i_Im_original,3);
    
    if rgb_check ==1
        v_rgbc = 0; %the image is greyscale
    else
        v_rgbc = 1; %the image is rgb
    end
    
    %cropping check
    if cropping == 1
        a_i_image_cropped = imcrop(a_i_Im_original, [crop_x1 crop_y1 crop_x2 crop_y2]);
        a_i_image_to_be_processed = a_i_image_cropped; %this prevents black image on the final figure if cropping occurs.
        
    else
        a_i_image_to_be_processed = a_i_Im_original;
    end
    
    if v_rgbc == 1
        a_i_image_to_be_processed = rgb2gray(a_i_image_to_be_processed);
    end
    
    
%     testing_input = true;
%     if testing_input == 1
%         figure
%         imshow(a_i_image_to_be_processed)
%         return
%     end
    %% *Image Segmentation*
    % The original image is segmented to remove particles from the background. The
    % result is a binary image of particle (white) and background (black).

    %%   Filter based on local standard deviation.
    % The particles are first separated from the background by separating the image
    % based on local standard deviation. The particles will typically have a relatively
    % high local standard deviation when compared to the background of the optical
    % micrograph. The neighborhood (NHOOD) and threshold are adjusted to separate
    % the particles from the background as much as possible.
    
    %The particles are first segmented from the background by thresholding the
    %particles based on the local standard deviation of the particles.

    a_i_stdfilt = stdfilt(double(a_i_image_to_be_processed), ones(stdfilt_NHOOD)); %find the local stdev of the original image
    a_i_rough_filt = a_i_stdfilt > stdfilt_threshold; %threshold the local stdev. Creates a mask.
    a_i_rough = a_i_image_to_be_processed .* uint8(a_i_rough_filt); %Mask the original image
    
%     show_stdfilt_image = true;
%     if show_stdfilt_image ==1
%         figure
%         imshow(a_i_rough);
%         return
%     end

    %% A Gaussian filter is applied to blur the particles.
    a_i_image_gauss = imgaussfilt(a_i_rough,gaussian_sigma); %2D Gaussian image smoothing.
    
%     show_gaussian_image = true;
%     if show_gaussian_image == 1
%         figure
%         imshow(a_i_image_gauss);
%         return
%     end
    %%   Adjust the image intensities.
    % The particles are adjusted from their intensities in the Gaussian-filtered
    % image (LOW_IN & HIGH_IN) to a new range (LOW_OUT & HIGH_OUT).  This is done
    % to prepare the image for a  threshold, which is applied to segment the more
    % emphasized particles from the background further.
    
    %The image intensities are adjusted from intensity IN range from the
    %Gaussian filtered image, to an adjusted intensity image OUT.  This is done
    %to further emphasize the particles from the background.
    
    v_im_adj_in = [image_adjust_LOW_IN image_adjust_HIGH_IN]/255; %intensity values to be remapped
    v_im_adj_out= [image_adjust_LOW_OUT image_adjust_HIGH_OUT]/255; %intensity values remapped to
    
    a_i_image_adj= imadjust(a_i_image_gauss, v_im_adj_in, v_im_adj_out); %remap particles intensities to increase contrast.
%     
%     show_adjusted_image = true;
%     if show_adjusted_image == 1
%         figure
%         imshow(a_i_image_adj);
%         return
%     end
    %%   Threshold the image intensity.
    % A threshold is then applied to the newly adjusted image to create a binary
    % mask of particles features vs. background.  The threshold value can be determined
    % by inspecting the image intensities of the particles in the previous adjusting
    % step.
    
    % A threshold is applied to create a mask of the original image as
    % particles vs. non particle features.
    
    a_i_image_bw = a_i_image_adj > image_intensity_threshold; %color value greater than x to define particle as white. 150 thru 200; 170
%     
%     show_intensity_image = true;
%     if show_intensity_image ==1
%         figure
%         imshow(a_i_image_bw);
%         return
%     end
    
    %Fill in image holes
%     a_i_image_bw = imfill(a_i_image_bw, 'holes'); %enclose particle space if there are holes in the mask.
    %                 imshow(a_i_image_bw);
    %                 return
    
    
    %%   Filter particles based on size.
    % Non-particle features are further removed from the image by filtering based
    % on size.  This removes both large and small leftover regions from the background.
    % The smallest particle area (low) and largest particle area (high) define the
    % range of regions to be kept.  All features out of range are removed to the background.
    
    %The particle mask is further filtered based on region size.  Small and
    %large non particle features are removed.
    
    a_i_image_filtered = bwareafilt(a_i_image_bw, [area_filter_low area_filter_high]); %remove non-particles
%     
%     show_sizefiltered_image = true;
%     if show_sizefiltered_image ==1
%         figure
%         imshow(a_i_image_filtered);
%         return
%     end
%     
    
    %%   Filter particles based on eccentricity.
    % The image is further segmented based on the eccentricity of the particles.
    % Some features have an area that falls within the pre-defined range, but still
    % need to be removed. The eccentricity of each region is calculated, and a threshold
    % is applied to remove non-particle features based on eccentricity.
    
    %The remaining regions are labeled as connected, and the eccentricity of
    %each region is calculated.  A threshold is then applied to segment the
    %image further based on eccentricity, which removes non-particle features
    %that have an area that was within the acceptable range.
    
    [a_i_labeled_image, i_particles_count_eccentricity] = bwlabel(a_i_image_filtered); %determine connected regions.
    regionprops_eccen = regionprops(a_i_labeled_image, 'Perimeter','Eccentricity', 'PixelIdxList'); %calculate eccentricity
    
    % eccentricity_threshold = 0.2; %eccentricity threshold.
    
    % a_c_particles_to_remove = [c_regionprops_eccen.Eccentricity] < eccentricity_threshold ; % Apply threshold, may need to flip sign depending on the image.
    % a_i_image_eccen = a_i_image_filtered;
    % a_i_image_eccen( vertcat( c_regionprops_eccen(a_c_particles_to_remove).PixelIdxList ) ) = 0; % remove particles who are members of 'toremove'
    

particles_to_remove = [regionprops_eccen.Eccentricity] < eccentricity_threshold | [regionprops_eccen.Perimeter] > perimeter_threshold ; % Apply threshold, may need to flip sign depending on the image.
a_i_image_eccen = a_i_image_filtered;
a_i_image_eccen( vertcat( regionprops_eccen(particles_to_remove).PixelIdxList ) ) = 0; % remove particles who are members of 'toremove'

%     show_eccentricity_image = false;
%     if show_eccentricity_image ==1
%         figure
%         imshow(a_i_image_eccen);
%         return
%     end
    
    
    % figure; imshow(a_i_image_eccen);
    % figure; imshow(a_i_image_to_be_processed);
    % return
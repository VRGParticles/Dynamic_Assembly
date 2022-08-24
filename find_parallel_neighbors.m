function       [para_accepted_particles, para_accepted_sum, para_kept] = find_parallel_neighbors(particle_distances, angle_differences, param);


% ; NAME:
% ; find_parallel_neighbors
% ; PURPOSE:
% ; Find which particles are parallel to each other (major axis)
% ; background.
% ; CATEGORY:
% ; Particle Classification 
% ; CALLING SEQUENCE:
% ; [para_accepted_particles, para_accepted_sum, para_kept] = find_parallel_neighbors(particle_distances, angle_differences, param);
% ; INPUTS:
% ; particle_distances: the center to center distance between particles.
% ; angle_differences: the cdifference in angle between the major axes of
% the particles
%   param:  a structure containing parameters for the various parts of the
%   image editing process: cropping, filtering, etc. Each is defined in the comments of the parameter
%   file.
% ; OUTPUTS:
% ; para_accepted_particles: logical array showing which particle pairs are
% counted
% ; para_accepted_sum: summation of how many neighbors a particle is
% parallel to.
% ; para_kept: logiacl array showing if the particle is counted as a
% parallel member


 %for the chains, we consider parallel angles.
    para_angle1_low = param.para_angle1_low; %these are given as examples, will not be true if absolute value of the angles are considered.
    para_angle1_high = param.para_angle1_high;
    
    %second set of angles for set 2.
    para_angle2_low = param.para_angle2_low;
    para_angle2_high = param.para_angle2_high;
    
    particle_distance_threshold= param.particle_distance_threshold;
    
    
    accepted_distances = (0 <= particle_distances & particle_distances <= particle_distance_threshold); %must be within threshold for distance between centers
    
    para_accepted_angles = (para_angle1_low <= angle_differences ...
        & angle_differences <= para_angle1_high) ...
        | (para_angle2_low <= angle_differences ...
        & angle_differences <= para_angle2_high);
    
    para_accepted_particles = (accepted_distances == para_accepted_angles) & (accepted_distances ==1); %determine if particles can be accepted by para if they meet threshold
    
    para_accepted_sum=sum(para_accepted_particles,2); % move to one column, sums how many times the particles meet the conditions.
    
    para_kept = para_accepted_sum >= 1;%logical classification if the particle is classified (kept) as para or not.
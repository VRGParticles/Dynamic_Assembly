function angle_differences = find_angle_differences(particle_distances, number_of_particles, regionprops_final, param);

% ; find_angle_differences
% ; PURPOSE:
% ; find the angle differences between the major axes of the particles.
% ; CATEGORY:
% ; Particle Classification
% ; CALLING SEQUENCE:
% ; angle_differences = find_angle_differences(particle_distances, number_of_particles, regionprops_final, param)
% ; INPUTS:
% ; particle_distances: the center to center distance between particles.
% ; number_of_particles: the number of particles detected
% ; regionprops_final: the properties of each particle
% ;  param:  a structure containing parameters. Each is defined in the comments of the parameter
%   file.
% ; OUTPUTS:
% ; angle_differences: an array containing the angle differences between particles i and j that meet the defined distance threshold and pass other tests.

%Preallocate arrays
angle_differences = ones(number_of_particles)*-1000; %We are setting the preallocation as -1000 in each cell, as if we have any reasonable angle(in degrees), the 'blank' value could potentially be picked up in the acceptable range. Will be prevented

%A distance threshold is applied to determine which particles are neighboring.  Then, the orientation angle difference of each particle is calculated.
particle_distance_threshold = param.particle_distance_threshold; %Threshold to determine the radius about each particle center where neighbors' centers fall.

%restrictions to count a particle as a potential neighbor.
malmin = param.malmin; %minimum length a particle must be (in major axis), set to 0 if you don't want this
areamin = param.areamin; %minimum a_c_area a particle must have, set to 0 if you don't want this
perimin = param.perimin; %minimum perimeter a particle must have. set to 0 if you don't want this
e2c =  param.e2c; %cutoff distance of ends of each particle to the center of a neighbor.

max_chain_length = param.max_chain_length;

[particle1, particle2] = find(particle_distances < particle_distance_threshold);  %determine which particles pairs are currently are accepted neighbors
chain_particle_indexes = cell(length(particle1), max_chain_length);
distcent1 = zeros(number_of_particles, number_of_particles);
distcent2 = zeros(number_of_particles, number_of_particles);


%Begin to calculate angle differences
for particle = 1:length(particle1)
    chain_particle_indexes{particle,1}=particle1(particle); %read in the ith particle
    chain_particle_indexes{particle,2}=particle2(particle); %read in the jth particle
    
    %map out the lines' end points and start points on Major Axis.
    xMajor1 = regionprops_final(chain_particle_indexes{particle,1}).Centroid(1) + (((regionprops_final(chain_particle_indexes{particle,1}).MajorAxisLength)./2) * cosd(regionprops_final(chain_particle_indexes{particle,1}).Orientation)); %x1
    xMajor2 = regionprops_final(chain_particle_indexes{particle,1}).Centroid(1) - (((regionprops_final(chain_particle_indexes{particle,1}).MajorAxisLength)./2) * cosd(regionprops_final(chain_particle_indexes{particle,1}).Orientation)); %y1
    yMajor1 = regionprops_final(chain_particle_indexes{particle,1}).Centroid(2) - (((regionprops_final(chain_particle_indexes{particle,1}).MajorAxisLength)./2) * sind(regionprops_final(chain_particle_indexes{particle,1}).Orientation)); %x2
    yMajor2 = regionprops_final(chain_particle_indexes{particle,1}).Centroid(2) + (((regionprops_final(chain_particle_indexes{particle,1}).MajorAxisLength)./2) * sind(regionprops_final(chain_particle_indexes{particle,1}).Orientation)); %y2
    
    %Calculate end to center distances of neighboring particles.
    distcent1(chain_particle_indexes{particle,1},chain_particle_indexes{particle,2})=pdist2([xMajor1 yMajor1],[regionprops_final(chain_particle_indexes{particle,1}).Centroid(1) regionprops_final(chain_particle_indexes{particle,1}).Centroid(2)],'euclidean'); %Point (x1,y1) of particle 1 to the center of particle 2.
    distcent2(chain_particle_indexes{particle,1},chain_particle_indexes{particle,2})=pdist2([xMajor2 yMajor2],[regionprops_final(chain_particle_indexes{particle,2}).Centroid(1) regionprops_final(chain_particle_indexes{particle,2}).Centroid(2)],'euclidean'); %Point (x2,y2) of particle 1 to the center of particle 2.
    
    %Use restrictions to see if the particle should be counted as a
    %chain member.  To be discounted, the particle must meet ALL
    %restrictions for the MAL, perimeter, and area (i.e. too small in
    %those three fields).  After it passes those restrictions, it moves
    %on to the ends to center distances check with the particle
    %neighbor.  If they do not meet the e2c threshold, they are not counted.  This
    %prevents particles that are not in a chain and are instead
    %adjacent but ajar from being counted as chain members.
    
    if (regionprops_final(chain_particle_indexes{particle,1}).MajorAxisLength > malmin || regionprops_final(chain_particle_indexes{particle,1}).Area > areamin || regionprops_final(chain_particle_indexes{particle,1}).Perimeter > perimin) ...
            && (distcent1(chain_particle_indexes{particle,1},chain_particle_indexes{particle,2}) < e2c) && (distcent2(chain_particle_indexes{particle,1},chain_particle_indexes{particle,2}) < e2c)  %ensure the particle is not too small based on a_c_area and perimeter
        
        angle_differences(chain_particle_indexes{particle,1}, chain_particle_indexes{particle,2}) = abs(regionprops_final(chain_particle_indexes{particle,1}).Orientation - regionprops_final(chain_particle_indexes{particle,2}).Orientation); %difference in angles'
    end
end


for particle=1:number_of_particles
    angle_differences(particle, particle) = -10000; %remap the distance between the same particle (1 to 1, for example) from 0 to -1000.  This is to keep these particles from inflating the neighbor count, as '0' can potentially fall into an angle range.
end
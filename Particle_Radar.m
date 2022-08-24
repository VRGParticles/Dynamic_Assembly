%%Particle Radar- Colloidal Chain Counting in optical micrograph videos.

% This script is useful for quantifying order in optical micrographs of 2-D anisotropic colloidal particle assemblies.
% *Input:* *image* (file extension .png, .jpg, .tif, etc.) as the optical
% micrograph to be analyzed. The images should be in a stack of frames,
% named chronologically. In this script, the file names are:
% 'frame_cropped_#', where # is the numerical frame number.

%Each frame is read in and passed to the classifier, which locates the
%particles and counts the chains. The location (centroids) are collected
%and passed to the tracker (Trackpy). The tracking data is then read back
%to MATLAB, merged with the classification data, then analyzed.

%most parameters are in the 'Parameters_Particle_Radar.m' file.
%% 

%     Copyright (C) 2022 Veronica Grebe, Weck Group at New York University.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see http://www.gnu.org/licenses/

disp('Copyright (C) 2022 Veronica Grebe, Weck Group at New York University')
disp('This program comes with ABSOLUTELY NO WARRANTY.')
disp('This is free software released under GNU GPL-3.0,')
disp('and you are welcome to redistribute it under certain conditions;')
disp('see http://www.gnu.org/licenses/')
%% Setup
positionindex = 1; %keeping track of writeout
tracking_position_inputs = zeros(20*10, 3); %inputs that will go to the tracker
positions_per_frame = cell(1,1); %all of the positions kept.
classification_per_frame = cell(1,1); %all of the classifications kept.

% frame read in
howmanyframes = 450; %total number of frames in the stack that you want to at least run through the classifier portion.
chains_per_frame = cell(7, howmanyframes); %all data per all frames, see descriptions below as assigned

run('Parameters_Particle_Radar.m'); %run the initial parameters
cd Frames %change directory to wherever your image stack is

%% BEGIN ANALYSIS
%The image is read in, detected as greyscale/rgb, and cropped as determined
%by the user.

for videoFrame = 1:howmanyframes
    a_i_Im_original = imread(['frame_cropped_', num2str(videoFrame), '.png']);

    if param.cropping == 1
        image = imcrop(a_i_Im_original, [param.crop_x1 param.crop_y1 param.crop_x2 param.crop_y2]);
    end

    %Preprocess the image/video frame and convert it to binary
    a_i_image_eccen = image_preprocess(a_i_Im_original, param);

    %%  Locate and define the properties of the fully preprocessed image.
    [image_connected,number_of_particles] = bwlabel(a_i_image_eccen); %show connectivity of final preprocessed image

    regionprops_final = regionprops(image_connected,'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'PixelIdxList', 'Perimeter'); %obtain properties for the connected regions

    cns = zeros(number_of_particles,2); %prealllocate centroids. column 1 is x, 2 is y.

    %write out the centroids
    for particle = 1:number_of_particles
        cns(particle,:) = regionprops_final(particle).Centroid;
    end
    %% Prepare tracker input  particle positions and frames

    %put into order for tracker. column 1 is x centroid, 2 is y centroid,
    %and 3 is frame number.
    tracking_position_inputs(positionindex:(positionindex + number_of_particles -1),1) = cns(:,1); % x values
    tracking_position_inputs(positionindex:(positionindex + number_of_particles -1),2) = cns(:,2); % yvalues
    tracking_position_inputs(positionindex:(positionindex + number_of_particles -1),3) = videoFrame; %Frame number


    positionindex = positionindex + number_of_particles; %increase current index to write out to in the tracker position inputs

    %% Strip out the locations in each frame for later searching
    %To ensure matching during the inner join
    %later, we have to 'fix' the values by truncating.
    positions_per_frame{1,videoFrame} = [fix(cns(:,1)*10^4)/10^4, fix(cns(:,2)*10^4)/10^4];


    %% *Classification Rules* Distance and Orientation Calculations
    % The Euclidean distance (d) and difference in orientation angle (?) between
    % particles are calculated to setup classification rules.  Particles are counted
    % as non-chain members based of restriction criteria set by the user.

    % center-to-center distance calculation
    particle_distances = squareform(pdist(cell2mat({regionprops_final(:).Centroid}'))); %computes Euclidean distance from particle i centroid to particle j centroid (also known as d).


    %%   Find particles within a distance threshold (D) and within orientation threshold.
    % The Euclidean distance between particle centroids is calculated.  To determine
    % which particles are neighboring, a threshold of distance is applied, which projects
    % out from particle centroids as a radius.  Two particles whose centroids fall
    % within the threshold are considered neighbors.
    angle_differences = find_angle_differences(particle_distances, number_of_particles, regionprops_final, param);

    %% Angle relationship: Parallel.
    % The angle relationship that is used to classify particles as in a chain is
    % parallel (abbr. para).  The range is extended from 0 (and 180) degrees to give
    % tolerance in the particle assemblies. We used 30 degrees as the tolerance.
    % Two ranges are accepted for the parallel range.  These numbers should be adjusted - along with the distance threshold
    % - to ensure that the amount of neighbors that fall within the angle relationship
    % range are as expected.

    %% Find which particles meet the angle AND distance requirement.
    [para_accepted_particles, para_accepted_sum, para_kept] = find_parallel_neighbors(particle_distances, angle_differences, param);
    %para_accepted_particles is a connectivity matrix showing where
    %particles [row,column] have an accepted angle relationship (1). Rejected are (0)
    %para_accepted_sum shows a total count of how many particles each particle has
    %the relationship with
    %para_kept is a logical array showing if the particle [row] has hit the
    %angle requirement at least once.

    %If you want to visualize the particles, uncomment below. Note the
    %return statement at the end.

    % % %    plot_single_frame_angle_counts(a_i_Im_original, param, number_of_particles, regionprops_final, para_accepted_sum); return


    %% Count the chains
    total_chain_count = zeros(1,param.max_chain_length); %Where the results of the total chain count will go.

    %addressing when neighbors >2 at some point due to aggregation or
    %branching. this section is more whimsical than others, have to check
    %results or choose to run this function or not.

    %Chain correction
    para_accepted_particles = branched_chains_fix(para_accepted_particles, regionprops_final, number_of_particles, param);

    %     figure; plot(digraph(para_accepted_particles)); % plot directed graph
    G = digraph(para_accepted_particles);
    bins = conncomp(G, 'Type','weak'); %weak means one direction. The content is the chain id #. if you count the # of times the particle id number shows up that is the chain length. max number is # of chains.
    bins = bins';
    number_of_chains = max(bins); %total number of chains
    chain_counts_accum = zeros(30,1); %accumulating term
    chain_counts =zeros(number_of_particles, 1); % which particle is what number

    %Count the number of each chain type
    for chain_number = 1:number_of_chains
        count = find(bins(:) == chain_number); %which particles are part of chain
        chain_counts_accum(length(count)) = chain_counts_accum(length(count)) + 1; %add one to the count
        for particle = 1:length(count)
            chain_counts(count(particle)) = length(count);
        end
    end

    %preallocation
    chain_list_all = cell(number_of_particles,1); %the particles in each chain assigned per partilce
    chain_list = cell(number_of_chains,1); %a list of all unique chains

    %build up the particle #'s in each chain
    for chain = 1:number_of_chains
        chain_list{chain,1} = find(bins == chain); %chain number is the row #, and the particles in it are in the list by local id number.
    end


    for particle = 1:number_of_particles
        chain_list_all(particle,1) = chain_list(bins(particle,1)); %chain_list but written to particle instead of chain
    end

    %Collect all of the classification information
    chains_per_frame{1,videoFrame} = chain_list; %lists of chains
    chains_per_frame{2,videoFrame} = bins; %which chain number the particle is in
    chains_per_frame{3,videoFrame} = chain_list_all; %every particle full chain list, chain_list but written to particle instead of chain
    chains_per_frame{4,videoFrame} = [1:number_of_particles]'; %number of particles
    classification_per_frame{1,videoFrame} = chain_counts; %write out the length of the chain
    positions_per_frame{1,videoFrame}(:,3) = chain_counts;

    %kinetics calculation in the frame, for number-averaged degree of
    %polymerization.

    numerator = 0;
    denominator = length(chains_per_frame{1,videoFrame}); %how many total chains there are
    for chains = 1:length(chains_per_frame{1,videoFrame}) %all chains
        numerator = numerator + length(cell2mat(chains_per_frame{1,videoFrame}(chains,1))); %add in all particles in chains to numerator
    end

    %assign kinetic dataj
    chains_per_frame{5,videoFrame} = numerator;
    chains_per_frame{6,videoFrame} = denominator;
    chains_per_frame{7,videoFrame} = numerator/denominator;

    %% Plot Chain count of single frame. 
    %uncomment below to plot the chain count of a single frame. Note the
    %return statement at the end.
% % %     plot_single_frame(a_i_Im_original, param, number_of_particles, regionprops_final, chain_counts, chain_counts_accum)
% % %     return
end

%%  Pass particle positions to the tracker
%tr is the result
trT = array2table(tracking_position_inputs, 'VariableNames', {'x','y', 'frame'});
writetable(trT,'Tracker_Input.csv')
disp('Read in Tracker_Input.csv into Trackpy now to track particles.')
disp('Strike any key to continue and read in tracked data.')
pause
%% duration
[tr tr_sorted results_sorted_particles] = process_tracking_data(howmanyframes, positions_per_frame, chains_per_frame);

disp('Data merged. For ideas for analysis, see analysis section of the program.')
return

%The merged results are in results_sorted_particles as a labeled table. 
%In the unlabeled array tr_sorted, the columns are as follows:
%1 - x position
%2- y position
%3 - frame #
%4 - particle #
%5 - assigned chain count
%6 - local ID
%7 - local chain number

%% NOTE:
%Below are the different types of analysis that can be done on the tracker
%data. 
%% tracker write off based on classification color
%create colors as outlines based on particle ID numbers.

colors = load('colors.mat'); %RNG colors
colors = colors.colors;

colors = vertcat(colors, colors); 
colors = vertcat(colors, colors); 
colors = vertcat(colors, colors); 
colors = vertcat(colors, colors); 
colors = vertcat(colors, colors); 
colors = vertcat(colors, colors); 
colors = vertcat(colors, colors); 

%bar colors to plot bar graphs, can change
bar_colors = {[0 0 0] [1 0 0], [1 1 0], [0 1 0], [0 1 1],...
    [0 0 1], [0.54 0.17 0.89] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] ...
    [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] ...
    [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] ...
    [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] ...
    [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] [1 0 1] ...
    }; % set up a color palette for the legend labels later

bar_colors = horzcat(bar_colors, bar_colors);
%% Plot the tracks of the particles on each frame
% label border and tracks are basd on particle ID, and the label itself is
% colored from the class of the particle.

for framenumber = 1:howmanyframes %customize
    
    image = imread(['frame_cropped_', num2str(framenumber), '.png']); %read in the frame
    if param.cropping == 1
        image = imcrop(image, [param.crop_x1 param.crop_y1 param.crop_x2 param.crop_y2]);
        f = figure('visible','off'); imshow(image);
    else
        f = figure('visible','off'); imshow(image);
    end

    hold on

    disp(['Now processing frame number ' num2str(framenumber)])

    frame_rows = zeros(1,1); %preallocate/reset

    frame_rows = find(tr_sorted(:,3) == framenumber); %where tr_sorted rows contains the current frame number

    if framenumber == 1 %special case where we just want a static plot

        for particle = 1:length(frame_rows)
            scatter(tr_sorted(frame_rows(particle),1), tr_sorted(frame_rows(particle),2), 11, 'MarkerFaceColor', bar_colors{1,tr_sorted(frame_rows(particle),5)}, ...
                'MarkerEdgeColor', colors{tr_sorted(frame_rows(particle),4),1}, 'MarkerEdgeAlpha',1, 'LineWidth',1)
        end
        export_fig(['plots_big',num2str(framenumber) '.tif'], '-native', '-tif');

    else

        for particle = 1:length(frame_rows) %go through all current frame rows currently in framerate

            rows = zeros(2,1); %reset
            all_id_rows = zeros(1,1); %reset
            all_frames_rows = zeros(1,1);
            last_frame = [];
            rows(2,1) = frame_rows(particle);%the higher frame rate row id is in row 2.
            all_id_rows = find(tr_sorted(:,4) == tr_sorted(frame_rows(particle),4));
            all_frames_rows = find(tr_sorted(:,3) == (framenumber - 1)); %all frames before it
            last_frame = intersect(all_id_rows, all_frames_rows);

            if isempty(last_frame) == 0 %last frame found
                rows(1,1) = last_frame;
            end

            %If only one id exists, the current frame, so just show one point instead
            %of the connected line.
            if rows(1,1) == 0 %the particle id does not exist in the previous frame
                scatter(tr_sorted(frame_rows(particle),1), tr_sorted(frame_rows(particle),2), 11, 'MarkerFaceColor', bar_colors{1,tr_sorted(frame_rows(particle),5)}, ...
                    'MarkerEdgeColor', colors{tr_sorted(frame_rows(particle),4),1}, 'MarkerEdgeAlpha',.5, 'LineWidth',1)
            end

            if rows(1,1) > 0 %current particle also existed in the last frame.

                for rownumber = 1:length(rows)
                    if rownumber < length(rows) %the previous frame
                        plot(tr_sorted(rows,1), tr_sorted(rows,2), '-', 'LineWidth',1,'MarkerSize', 2, 'color', colors{tr_sorted(frame_rows(particle),4),1});

                        scatter(tr_sorted(rows(1,1),1), tr_sorted(rows(1,1),2), 11, 'MarkerFaceColor', bar_colors{1,tr_sorted(rows(1,1),5)}, ...
                            'MarkerEdgeColor', colors{tr_sorted(frame_rows(particle),4),1}, 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.6, 'LineWidth',1)
                    end
                    if rownumber == length(rows) %current frame
                        scatter(tr_sorted(rows(2,1),1), tr_sorted(rows(2,1),2), 11, 'MarkerFaceColor', bar_colors{1,tr_sorted(frame_rows(particle),5)}, ...
                            'MarkerEdgeColor', colors{tr_sorted(frame_rows(particle),4),1}, 'MarkerEdgeAlpha',1,'LineWidth',1)
                    end
                end

            end
        end

        export_fig(['plots_big',num2str(framenumber) '.tif'], '-native', '-tif');

    end
end

%% Plot the duration for each particle, meaning for how many frames do they persist

for particle = 1:max(tr(:,4))
    dura(particle,1) = length(find(tr(:,4) == particle));
end
figure; histogram(dura)
xlabel('Amount of Frames', 'Fontsize', 10)
ylabel('Amount of Particles', 'Fontsize', 10)
%% Plot the Duration on a single image
for framenumber = 200 %customize
    %
    duration = 300; %minimum amount of frames the particle must have been detected in
    hits = 0;
    %bc we have 10 in testing
    image = imread(['frame_cropped_', num2str(framenumber), '.png']); %read in the frame
    if param.cropping == 1
        image = imcrop(image, [param.crop_x1 param.crop_y1 param.crop_x2 param.crop_y2]);
        f = figure; imshow(image);
    else
        f = figure; imshow(image);
    end

    hold on

    frame_rows = zeros(1,1); %preallocate/reset

    frame_rows = find(tr_sorted(:,3) == framenumber); %where tr_sorted rows contains the current frame number


    for particle = 1:length(frame_rows) %go through all current frame rows currently in framerate
        if dura(tr_sorted(frame_rows(particle),4)) > duration
            hits = hits + 1; %how many particles with the given duration are in the frame

            rows = zeros(2,1); %reset
            all_id_rows = zeros(1,1); %reset
            all_frames_rows = zeros(1,1);
            last_frame = [];
            rows(2,1) = frame_rows(particle);%the higher frame rate row id is in row 2.
            all_id_rows = find(tr_sorted(:,4) == tr_sorted(frame_rows(particle),4));
            all_frames_rows = find(tr_sorted(:,3) == (framenumber - 1)); %all frames before it
            last_frame = intersect(all_id_rows, all_frames_rows);

            %add text label
            text(tr_sorted(frame_rows(particle),1) -6, tr_sorted(frame_rows(particle),2),num2str(dura(tr_sorted(frame_rows(particle),4))), ...
                'color', [.19 .98 .19], 'FontSize', 8 )


        end
    end


    saveas(f, ['plot_labeled_long_', num2str(framenumber), '.png']);
    disp([ 'total number of particles that meet duration in this frame: ' num2str(hits) ] )
end


%% Plot frame-by-frame duration changes.

tr_sorted_for_dura = table2array(results_sorted_particles(:,1:7)); %make a new result table of only duration met particles
tr_dura = []; %preallocate
row_new = 1; %index of the new array
duration = 300;
for rows = 1:length(tr_sorted) %all rows
    if dura(tr_sorted_for_dura(rows,4)) >= duration %particle in rows has a big enough dura
        tr_dura(row_new,:) = tr_sorted_for_dura(rows,:); %import into tr_new all columns of tr_sorted
        row_new = row_new + 1;
    end
end

dura_particles = unique(tr_dura(:,4)); %particle ids that meet dura. These are the particles we will analyze in the tracker for now on.
%% plot one particle over time, sample!

figure; plot( tr_dura(1:445 ,3) , tr_dura(1:445, 5), ':O', 'color','black', 'LineWidth', 1.5, ...
    'MarkerEdgeColor', [0 0 139/255], 'MarkerFaceColor', [173/255 216/255 230/255])
set(gca,...
    'FontSize', 30, ...
 'XAxisLocation','bottom',...
 'YAxisLocation','left',...
 'XColor','k','YColor','k',...
 'LineWidth',2,...
 'TickDir','out', ...
 'box', 'off');

xlabel('Frame Number', 'Fontsize', 30)
ylabel('Chain Length', 'Fontsize', 30)
ylim([0 10])
disp('Done')




%% Count the amount of each chain lenght, then convert to percent

sample_total = unique(tr_dura(1:445,5)); %find the unique chain ids it hit
for row = 1:length(sample_total)
    sample_count(sample_total(row,1),1) = numel(find(tr_dura(1:445,5) == sample_total(row,1))); %count them
end

for row = 1:length(sample_count)
    sample_count(row, 2) = (sample_count(row,1)/ sum(sample_count(:,1))) * 100; %convert to percent
end

bars = [0 0 0;1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;0.54 0.17 0.89;1 0 1; 1 0 1];

figure;
b = bar(sample_count(:,2), 'FaceColor', 'flat');
b.CData = bars;
ylabel('Percent of Frames', 'Fontsize', 30)
xlabel('Chain Length', 'Fontsize',30)
set(gca,'TickDir', 'out');
set(gca,'FontSize', 30)
set(gca,'LineWidth', 2, 'box','off')


%% How do the particle chains change
% row # is the PREVIOUS frame chain length (frame >=2, looks at frame -1)
%MAJOR ASSUMPTION: CHANGES BETWEEN FRAMES ARE ALL THE SAME WEIGHT, EVEN IF
%THERE ARE 20 FRAMES IN BETWEEN, AS LONG AS IT MEETS DURATION.
%columns 1:7 -3+ -2 -1 0 1 2 3.

dura_results = zeros(20,7);
particle_id = 1; %the row index of dura_particles, preallocate
row = 1; %start at the second row
while  row < length(tr_dura)
    row = row + 1; %increment the current row being searched.

    if tr_dura(row,4) == tr_dura(row -1 ,4) %the particle id is identical
        change = tr_dura(row,5) - tr_dura(row -1 ,5); %change in the chain

        if change <= -3 %the chain breaks alot
            dura_results(tr_dura(row -1 ,5),1) = dura_results(tr_dura(row -1 ,5),1) + 1;
        elseif change == -2
            dura_results(tr_dura(row -1 ,5),2) = dura_results(tr_dura(row -1 ,5),2) + 1;
        elseif change == -1
            dura_results(tr_dura(row -1 ,5),3) = dura_results(tr_dura(row -1 ,5),3) + 1;
        elseif change == 0
            dura_results(tr_dura(row -1 ,5),4) = dura_results(tr_dura(row -1 ,5),4) + 1;
        elseif change == 1
            dura_results(tr_dura(row -1 ,5),5) = dura_results(tr_dura(row -1 ,5),5) + 1;
        elseif change == 2
            dura_results(tr_dura(row -1 ,5),6) = dura_results(tr_dura(row -1 ,5),6) + 1;
        elseif change >= 3
            dura_results(tr_dura(row -1 ,5),7) = dura_results(tr_dura(row -1 ,5),7) + 1;
        end

    elseif tr_dura(row,4) ~= tr_dura(row -1 ,4) %particle id is changing
        particle_id = particle_id + 1 ;%increment the row for the particle id
    end
end

%% Convert to percent
dura_results_sums = zeros(18,7);
for row=1:length(dura_results_sums)
    row_sum = sum(dura_results(row,:));
    for column = 1:7
        dura_results_sums(row,column)=(dura_results(row,column)/row_sum) * 100;
    end
end

%% Make histogram of all
figure;
hold on;
for row = 1:length(dura_results_sums)
    bar(row, dura_results_sums(row,:), 'stacked')
end

xlabel('Chain Length in Previous Frame','Fontsize', 20)
ylabel('Percent of Change', 'Fontsize', 20)
ylim([0 100]);
ylim([0 100]);
set(gca,'TickDir', 'out',  'XColor','k','YColor','k',...
 'LineWidth',2);
set(gca,'FontSize', 20)

h = zeros(7, 1);
h(1) = plot(NaN,NaN,'ks', 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
h(2) = plot(NaN,NaN,'ks', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
h(3) = plot(NaN,NaN,'ks', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
h(4) = plot(NaN,NaN,'ks', 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
h(5) = plot(NaN,NaN,'ks', 'MarkerFaceColor', [0.9290 0.6940 0.1250]);
h(6) = plot(NaN,NaN,'ks', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
h(7) = plot(NaN,NaN,'ks', 'MarkerFaceColor', [0 0.4470 0.7410]);


legend(h, '+3 or more','+2','+1', 'No change', '-1', '-2', '-3 or less');

% legend(h, '-3 or less','-2','-1', 'no change', '+1', '+2', '+3 or more');

%% Kinetics plot
figure; scatter(linspace(1,90,450),[chains_per_frame{7,:}], 'MarkerFaceColor', 'blue')


set(gca,...
    'FontSize', 28, 'box', 'off', ...
 'XColor','k','YColor','k',...
 'LineWidth',2,...
 'TickDir','out');

xlabel('Time (Seconds)')
ylabel('Number-Averaged Degree of Polymerization')


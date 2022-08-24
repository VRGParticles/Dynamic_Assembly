function [tr tr_sorted results_sorted_particles] = process_tracking_data(howmanyframes, positions_per_frame, chains_per_frame)

% ; NAME:
% ; plot_single_frame
% ; PURPOSE:
% ; Plot the chain count in a single frame

% ; CALLING SEQUENCE:
% ; [tr tr_sorted results_sorted_particles] = process_tracking_data(howmanyframes, positions_per_frame, chains_per_frame)
% ; INPUTS:
% ; howmanyframes : number of frames analyzed
% ; positions_per_frame : structure that has all x and y values for every
% frame
% ; chains_per_frame : structure that has all chain data for all frames, as
% well as kinetic data.
% ; OUTPUTS:
% ; tr: unlabeled table of merged data
% ; tr_sorted: tracking results sorted by frame #
% ; results_sorted_particle : labeled table of merged data, sorted by
% particle number

tr = csvread('Tracked_Particles_Input.csv'); %edited, it is x y frame particle

for row = 1:length(tr)
    tr(row,4) = tr(row,4) + 1; %stop particle indexing at 0
end

%%
for row = 1:length(tr)
    tr(row, 5) = row;
end

max_test = max(tr(:,4));

%tr is currently: x,y,framenumber,particleid,row#
%% correct the rounding error 4 dp
for row = 1:length(tr)
    tr(row,1) = fix(tr(row,1)*10^4)/10^4; %stop particle indexing at 0
    tr(row,2) = fix(tr(row,2)*10^4)/10^4; %stop particle indexing at 0

end
%% join tables

%positions per frame rounding error fix
for iterations = 1:howmanyframes
    clear right_rows Tright Tleft Tinner_table Tinner
    %right side
    right_rows = find(tr(:,3) == iterations);
    Tright = table([tr(right_rows(1,1):right_rows(length(right_rows),1),1),tr(right_rows(1,1):right_rows(length(right_rows),1),2)], ...
        tr(right_rows(1,1):right_rows(length(right_rows),1),5),...
        'VariableNames',{'pos' 'row number'}); %x and y positions and row numbers. assuming they will never be equal.

    %left side
    Tleft = table(positions_per_frame{1,iterations}(:,3), [positions_per_frame{1,iterations}(:,1),positions_per_frame{1,iterations}(:,2)],...
        chains_per_frame{4,iterations}(:,1), chains_per_frame{2,iterations}(:,1), chains_per_frame{3,iterations}(:,1), 'VariableNames', ...
        {'chain_length' 'pos' 'localid' 'chain_number' 'chain_list'}); %chain length assigned per position. Two Placeholders for row and ID


    Tinner_table = innerjoin(Tleft,Tright); %inner join of the two

    Tinner = table2cell(Tinner_table); %cell


    for row = 1:length(Tinner)
        tr(Tinner{row,6},5) = Tinner{row,1}; % import the chain_length into the tracker table
        tr(Tinner{row,6},6) = Tinner{row,3}; %local id
        tr(Tinner{row,6},7) = Tinner{row,4}; %local chain number
    end
end

%convert the results table to a table
tr2 = num2cell(tr); %convert to cell to accept lists

%read in the chain list
for row = 1:length(tr2)
    tr2(row,8) = chains_per_frame{3,tr2{row,3}}(tr2{row,6},1); %chain list
end

results = cell2table(tr2,  'VariableNames', {'XPosition' 'YPosition' 'Frame' 'ParticleID' 'ChainCount' 'LocalID' 'LocalChainID' 'ChainList'});
results_sorted = sortrows(results,3); %sort on frame #
results_sorted_particles = sortrows(results,4); %sort on particle ID
tr_sorted = table2array(results_sorted(:,1:7)); %array, sorted on frame #
%results table is:
%1 - x position
%2- y position
%3 - frame #
%4 - particle #
%5 - assigned chain count
%6 - local ID
%7 - local chain number
%8 - chain list

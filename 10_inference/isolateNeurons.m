clear; close all; clc
% This program recovers neurons on which to infer from HDF5 file.


%% Load HDF5:

h5path = '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2019-03-26/Run 10/Analysis/HDF5/2019-03-26(Run10).h5';


%% Getting the neurons:

% Parameters:
keep = [97, 98, 132, 201, 223];
labels = h5read(h5path, '/Data/Brain/Labels');
coord = h5read(h5path, '/Data/Brain/ZBrainCoordinates');
% Coordinates parameters:
inferior_y = 0.63;
mid_x = 0.25;
inferior_x = 0.06;
% Final vector:
keepn = ones(size(labels, 1), 1);
% Algorithm:
for i = 1:length(keepn)
    % Checking labels:
    label_ok = find(labels(i, :) == 1)';
    if sum(sum(label_ok == keep)) == 0
        keepn(i) = 0;
        continue
    elseif coord(i, 2) < inferior_y && abs(coord(i, 1) - mid_x) < inferior_x
        keepn(i) = 0;
        continue
    end
end
% Lowering vector size:
keepn = find(keepn == 1);


%% Plotting (optional):

figure
scatter3(coord(keepn, 1), coord(keepn, 2), coord(keepn, 3))
axis equal
        

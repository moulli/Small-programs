clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte'))
% Matlab sheet to get vector labels that returns labels for each ZBrain
% Coordinates from HDF5 file, based on zlabels005.


%% Paths to files:

orientation = 'RAS';
h5_path = '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2019-03-26/Run 09/Analysis/HDF5/2019-03-26(Run09).h5';
zlabels_path = '/home/ljp/Science/Hippolyte/zlabels005.mat'; % name of object must be zlabels005


%% Function to return labels:

% Recover coordinates:
coord = h5read(h5_path, '/Data/Brain/Coordinates');
nneu = size(coord, 1);
% Create "false" ZBraingrid object:
method = 'False temp object';
zbrainsize = [0.496, 1.122, 0.276];
increment = 0.005;
gridsize = floor(zbrainsize ./ increment);
orientation = 'RAS';
zfalse = ZBraingrid(method, gridsize, orientation);
% Add coordinates:
stemp = struct;
stemp.name = 'temp'; stemp.path = h5_path; stemp.comment = 'temp';
stemp.orientation = orientation;
stemp.coordinates = coord;
stemp.correlation = ones(nneu, 1);
addDataset(zfalse, stemp);

% Load zlabels:
load(zlabels_path);
% Create labels vector:
regions = length(zlabels005);
labels = zeros(nneu, regions);

% Add labels for each brain region:
for i = 1:regions
    % Get regions of interest:
    ztemp = zlabels005(i);
    % Recover indexes of points:
    indtemp = ztemp.Zindex;
    % Compare to zfalse:
    booltemp = (zfalse.Zindex == indtemp');
    deltemp = any(booltemp, 2);
    % Takes associated neurons:
    neutemp = zfalse.Zneuron(deltemp, :);
    neutemp = neutemp(:);
    neutemp(neutemp == 0) = [];
    % Change right column of labels:
    labels(neutemp, i) = 1;
    fprintf('Iteration %.0f. \n', i);
end
    


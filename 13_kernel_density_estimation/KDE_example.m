clear; close all; clc


%% Load data

% path to data
path = {'2020-02-18(Run01).h5';
        '2020-02-21(Run01).h5';
        '2020-02-21(Run03).h5';
        '2020-02-27(Run02).h5'};
plen = length(path);
    
% pick regressor
keepneureg = 1;

% get coordinates
h5coords_path = '/Data/Brain/ZBrainCoordinates';
h5reg_path = strcat('/Data/Brain/Analysis/RegressionDFFNeuron/regDataTail/keepneureg', num2str(keepneureg));
coords = cell(plen, 1);
for p = 1:plen
    h5coo = h5read(path{p}, h5coords_path);
    h5reg = h5read(path{p}, h5reg_path);
    coords{p} = h5coo(h5reg, :);
end


%% Create KDE objects for each dataset

kde = cell(plen, 1);
for i = 1:plen
    kdei = KDE;
    kdei.add(coords{i});
    kde{i} = kdei;
end


%% Create KDE object including all datasets

kde = KDE;
for i = 1:plen
    kde.add(coords{i});
end


%% Access data in KDE

disp(kde.h) % bandwidth
disp(kde.dim) % dimension (should be 3)
disp(kde.len) % number of points
disp(kde.options.isovalues) % plot option: isovalues
disp(kde.options.increment) % plot option: increment in grid
disp(kde.options.gridsize) % plot option: size of grid


%% Modify bandwidth, isovalues, grid increment, gridsize (NOT NECESSARY)

kde.bandwidth(0.015);
kde.isovalues([0.01, 0.01, 0.6]); % as many isovalues as you want
kde.increment(0.02); % OR: kde.increment([0.01, 0.01, 0.005])
kde.gridsize([0.5, 1.2, 0.3]); % OR: kde.gridsize(1)


%% Plot

% get original bandwidth
kde.bandwidth(0.0128);

% you can plot using the properties saved in the KDE object
figure
kde.plot
axis equal
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')

% or you can specify parameters without changing them in KDE
figure
kde.plot('isovalues', [0.1, 0.2, 0.4], 'increment', 0.0175, 'gridsize', [0.496, 1.122, 0.276])
axis equal
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
% NB: if you want to modify the bandwidth, you have to do it internally.


%% VERY IMPORTANT

%  The probability density function can be evaluated at any point in space,
%  but the results is senseless. PDFs for continuous distributions are
%  supposed to return the probability of having points in a SUBSPACE of the
%  total space. For instance, if you evaluate a normal distribution of low
%  variance in 0, the results will be above 1. But if you evaluate de
%  probability to have a point between -0.001 and 0.001, then the result
%  will be lower that 1. And the sum over all space will be equal to 1.
%
%  Based on this fact that evaluating PDF on a particular point has no
%  probabilistic meaning, when using plot or contour, I made the choice to
%  divide all the grid by the maximum value of the PDF in the grid.
%  Therefore, ISOVALUES HAVE NO MATHEMATICAL GROUND. These are just
%  arbitrary values, with nothing to do with probabilities.
%
%  I added the function project in the KDE class, that returns these
%  projection, with no division whatsoever. These projections are the
%  actual values of the PDF.


%% Contour

% you can plot contour projected on a plan
figure
kde.contour(1) % you have to specify the direction of the projection
axis equal
xlabel('y-axis', 'Interpreter', 'latex')
ylabel('z-axis', 'Interpreter', 'latex')

% and you can specify parameters without changing them in KDE
figure
kde.contour(1, 'isovalues', [0.1, 0.2, 0.4], 'increment', 0.0175, 'gridsize', [0.496, 1.122, 0.276])
axis equal
xlabel('y-axis', 'Interpreter', 'latex')
ylabel('z-axis', 'Interpreter', 'latex')


%% Project

% you can obtain projections of PDF on a plan
project = kde.project(1);
figure
image(project, 'CDataMapping', 'scaled')

% and you can specify parameters without changing them in KDE
project = kde.project(1, 'increment', 0.0175, 'gridsize', [0.496, 1.122, 0.276]);
figure
image(project, 'CDataMapping', 'scaled')
% NB: as explained above, these are the actual PDF values, so isovalues is
% no longer a parameter.


%% An example of comparison
% Say you want to compare distributions of the different HDF5 to total
% distribution.

% each distribution
kde = cell(plen, 1);
for i = 1:plen
    kdei = KDE;
    kdei.add(coords{i});
    kde{i} = kdei;
end
figure
for i = 1:4
    subplot(2, 2, i)
    kde{i}.contour(2, 'isovalues', 0.05:0.05:0.9, 'increment', 0.015)
    axis equal
    xlabel('x-axis', 'Interpreter', 'latex')
    ylabel('z-axis', 'Interpreter', 'latex')
end

% total distribution
kde = KDE;
for i = 1:plen
    kde.add(coords{i});
end
figure
kde.contour(2, 'isovalues', 0.05:0.05:0.9, 'increment', 0.015)
axis equal
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('z-axis', 'Interpreter', 'latex')





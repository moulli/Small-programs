clear; close all; clc



%% Going to the right address:

cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-06-14/Run 03/Analysis/HDF5'
% cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-05-24/Run 08/Analysis/HDF5'



%% We can get the info of the file:

hippo = 'h5hippo.h5';



%% Based on this info, we can obtain the data inside:

coordinates = h5read(hippo, '/Data/RefCoordinates');
values = h5read(hippo, '/Data/Values');
stimulus = h5read(hippo, '/Data/Stimulus');
labels = h5read(hippo, '/Data/Labels');



%% Normalizing the data:

nstimulus = (stimulus-mean(stimulus)) ./ std(stimulus);
nvalues = (values-mean(values, 2)) ./ std(values, [], 2);



%% Computing covariance between neurons and stimulus:

covariance = sum(nvalues.*nstimulus, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimulus.^2, 2));



%% Separating positive and negative stimuli:

nstimabs = abs(nstimulus);
nstimpos = nstimulus .* (nstimulus > 0);
nstimneg = -nstimulus .* (nstimulus < 0);

covabs = sum(nvalues.*nstimabs, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimabs.^2, 2));
covpos = sum(nvalues.*nstimpos, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimpos.^2, 2));
covneg = sum(nvalues.*nstimneg, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimneg.^2, 2));



%% Linear approach between correlation to positive and negative stimuli:

% Normalizing data:
covp = (covpos - mean(covpos)) ./ std(covpos);
covn = (covneg - mean(covneg)) ./ std(covneg);
% Linear regression with Moore-Penrose inverse:
lcov = length(covp);
X = [ones(lcov, 1), covp];
MPinv = (X' * X) \ X' * covn;
covnguess = X * MPinv;
% Plotting linear aproximation:
figure
subplot(1, 2, 1)
plot(covp, covn, '.', covp, covnguess, 'k')
grid on
title('Correlation with positive and negative stimuli for each neurons with regression and confidence interval', 'Interpreter', 'latex')
xlabel('Neuron', 'Interpreter', 'latex')
ylabel('Correlation', 'Interpreter', 'latex')
hold on

% Plotting error, with intervalles de confiance:
intervalle1 = 3.63; % 1% confidence interval
intervalle2 = 5.15; % 0.2% confidence interval
intervalle3 = 7.19; % 0.05% confidence interval
error = sqrt(mean((covnguess - covn).^2));
fprintf('RMS error is %f \n', error);
higherror1 = find((covn > covnguess + intervalle1*error & covn <= covnguess + intervalle2*error) | (covn < covnguess - intervalle1*error & covn >= covnguess - intervalle2*error));
higherror2 = find((covn > covnguess + intervalle2*error & covn <= covnguess + intervalle3*error) | (covn < covnguess - intervalle2*error & covn >= covnguess - intervalle3*error));
higherror3 = find(covn > covnguess + intervalle3*error | covn < covnguess - intervalle3*error);
plot(covp(higherror1), covn(higherror1), '.', 'Color', [1, 0.9, 0.9])
plot(covp(higherror2), covn(higherror2), '.', 'Color', [1, 0.65, 0.65])
plot(covp(higherror3), covn(higherror3), '.', 'Color', [1, 0, 0])

% Plotting different neurons:
subplot(1, 2, 2)
scatter3(coordinates(higherror1, 1), coordinates(higherror1, 2), coordinates(higherror1, 3), [], [1, 0.9, 0.9], '.')
hold on
scatter3(coordinates(higherror2, 1), coordinates(higherror2, 2), coordinates(higherror2, 3), [], [1, 0.65, 0.65], '.')
scatter3(coordinates(higherror3, 1), coordinates(higherror3, 2), coordinates(higherror3, 3), [], [1, 0, 0], '.')
title('Neurons not in confidence interval', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
axis equal
view([0, 90])



%% Clustering these packs of neurons:

% addpath('~/Science/Hippolyte/clustering');
errkeep = [higherror3; higherror2];
virtcov = [ones(size(higherror3)); 0.5*ones(size(higherror2))];
coordkeep = coordinates(errkeep, :);
dmax = 0.03;
minpt = 10;
coreffect = 10;
linClus = dbscan3(coordkeep, dmax, minpt, virtcov, coreffect, 1);



%% Computing centroid for each of the found cluster:

numclus = length(unique(linClus)) - 1;
centroids = zeros(numclus, 3);
for i = 1:numclus
    centroids(i, :) = mean(coordkeep(linClus == i, :));
end
% Plotting centroid on the same figure as dbscan3 figure:
hold on
scatter3(centroids(:, 1), centroids(:, 2), centroids(:, 3), 200, 'filled', 'MarkerFaceAlpha', 7/8, 'MarkerFaceColor', [0, 0, 0])



%% Using getLabel function:

% Importing sparse matrix:
importfile('~/Science/Hippolyte/MaskDatabase.mat');
% Using function:
out = getLabel(centroids, coordinates, labels, 'ON', MaskDatabaseOutlines);





    


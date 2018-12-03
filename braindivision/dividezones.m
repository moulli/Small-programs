clear; close all; clc



%% Going to the right address:

% cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-06-14/Run 03/Analysis/HDF5'
cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-05-24/Run 08/Analysis/HDF5'



%% We can get the info of the file:

hippo = 'h5hippo.h5';



%% Based on this info, we can obtain the data inside:

coordinates = h5read(hippo, '/Data/RefCoordinates');
labels = h5read(hippo, '/Data/Labels');
values = h5read(hippo, '/Data/Values');
nneu = size(labels, 1);



%% Normalizing the data:

nvalues = (values-mean(values, 2)) ./ std(values, [], 2);



%% Importing MaskDatabase

addpath('~/Science/Hippolyte/Small-programs/inhibited_neurons');
importfile('~/Science/Hippolyte/MaskDatabase.mat');



%% Taking brain region for each neuron:

reg = [1, 94, 114, 260, 275];
labh = zeros(nneu, 1);
for i = 1:nneu
    labtemp = find(labels(i, :) ~= 0)';
    valtemp = sum(labtemp == reg);
    switch sum(valtemp)
        case 0
            labh(i) = 0;
        case 1
            labh(i) = reg(valtemp == 1);
        otherwise
            mlabtemp = mean(labtemp);
            [~, mintemp] = min((reg' - mlabtemp).^2);
            labh(i) = reg(mintemp);
    end
end

bzones = cell(5, 2);
mid = 0.25;
for i = 1:5
    for j = 1:2
        switch j
            case 1
                ktemp = ((labh == reg(i)) & (coordinates(:, 1) <= mid));
            case 2
                ktemp = ((labh == reg(i)) & (coordinates(:, 1) > mid));
        end
        bzones{i, j} = find(ktemp == 1);
    end
end



%% Constructing correlation matrix:

izone = 1;
jzone = 1;
valij = nvalues(bzones{izone, jzone}, :);
% lenbzones = size(bzones{izone, jzone}, 1);
% corij = zeros(lenbzones);
% tic
% for i = 1:lenbzones
%     corij(1:lenbzones-i, lenbzones+1-i) = sum(valij(1:lenbzones-i, :).*valij(lenbzones+1-i, :), 2) ./ sqrt(sum(valij(1:lenbzones-i, :).^2, 2).*sum(valij(lenbzones+1-i, :).^2, 2));
%     if mod(i, 10) == 0
%         fprintf('Iteration %.0f out of %.0f, in %.2f seconds \n', [i, lenbzones, toc]);
%     end
% end
covij = cov(valij');
covijb = covij - tril(covij);

outij = corClust(valij, 0.75);
maxij = max(outij(2, :));
valijnew = zeros(0, size(valij, 2));

        












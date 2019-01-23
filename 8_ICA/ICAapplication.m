clear; close all; clc



%% Going to the right address:

sinstim = 1;
switch sinstim
    case 0
        cd '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/2018-05-24/Run 08/Analysis/HDF5'
    case 1
        cd '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/2018-06-14/Run 03/Analysis/HDF5'
end



%% We can get the info of the file:

hippo = 'h5hippo.h5';



%% Based on this info, we can obtain the data inside:

coordinates = h5read(hippo, '/Data/RefCoordinates');
labels = h5read(hippo, '/Data/Labels');
values = h5read(hippo, '/Data/Values');
stimulus = h5read(hippo, '/Data/Stimulus');
nneu = size(labels, 1);



%% Normalizing the data:

nvalues = (values-mean(values, 2)) ./ std(values, [], 2);



%% Importing MaskDatabase

addpath(genpath('~/Science/Hippolyte/Small-programs'));
importfile('~/Science/Hippolyte/MaskDatabase.mat');



%% Taking brain region for each neuron:

% Taking labels for each neuron:
reglab = [275, 1, 94, 114, 260];
labelsnum = permute(labels .* (1:size(labels, 2)), [2, 3, 1]);
labelslab = permute(sum(labelsnum == reglab), [3, 2, 1]);
% Finding neurons misassigned:
problems = find(sum(labelslab, 2) ~= 1);
noproblems = find(sum(labelslab, 2) == 1);
lprob = length(problems);
labelslab(problems, :) = zeros(lprob, size(labelslab, 2));
% K-nearest neighbours to find new labels:
knn = 20;
coordstemp = coordinates(noproblems, :);
dist = sum((coordstemp - permute(coordinates(problems, :), [3, 2, 1])).^2, 2);
[~, mindist] = sort(squeeze(dist));
for i = 1:lprob
    labelstemp = sum(labelslab(noproblems(mindist(1:knn+1, i)), :));
    [~, maxlab] = max(labelstemp);
    labelslab(problems(i), maxlab) = 1;
end
% Gathering neurons in a cell:
bzones = cell(5, 2);
mid = 0.25;
for i = 1:5
    for j = 1:2
        switch j
            case 1
                ktemp = ((labelslab(:, i) == 1) & (coordinates(:, 1) <= mid));
            case 2
                ktemp = ((labelslab(:, i) == 1) & (coordinates(:, 1) > mid));
        end
        bzones{i, j} = find(ktemp == 1);
    end
end
% Plotting:
figure
hold on
for i = 1:5
    for j = 1:2
        scatter3(coordinates(bzones{i, j}, 1), coordinates(bzones{i, j}, 2), coordinates(bzones{i, j}, 3), [], rand(1, 3))
        axis equal
    end
end



%% ICA on the brain regions:

A = fast_ica(values(bzones{1, 1}, :), 100, 30);
% s = A \ values(bzones{1, 1}, :);
X = values(bzones{1, 1}, :);
[Am, Wm, sm] = fastICA(X, 20, round(size(X, 1)/20));
Xm = Am * sm;
if size(sm, 1) < size(X, 1)
    figure
    hold on
    for i = 1:size(sm, 1)
        plot(sm(i, :))
    end
end
% figure; subplot(2, 1, 1); hold on; for i = 1:500; plot(s(i, :)); end; subplot(2, 1, 2); hold on; for i = 1:500; plot(snew(i, :)); end

smnorm = (sm - mean(sm, 2)) ./ var(sm, [], 2);
figure
subplot(1, 2, 1)
plot(max(abs(smnorm), [], 2), '.')
subplot(1, 2, 2)
hist(max(abs(smnorm), [], 2), 30)
Amvar = Am .* var(sm, [], 2)';
Ammean = Am * mean(sm, 2);
Xmnew = Amvar * smnorm + Ammean;

figure
subplot(4, 4, [1:3, 5:7, 9:11])
image(log(abs(real(Amvar))), 'CDataMapping', 'scaled')
colorbar
subplot(4, 4, 13:15)
plot(sum(Amvar))
subplot(4, 4, 4:4:12)
image(Ammean, 'CDataMapping', 'scaled')
colorbar





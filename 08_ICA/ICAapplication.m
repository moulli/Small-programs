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



%% Neural network:

addpath('/home/ljp/Science/Guillaume/Spontaneous/Datasets')
% h5disp('20180706_Run04.h5')
coords = h5read('20180706_Run04.h5', '/Data/zbrain_coords');
coords(:, 2) = max(coords(:, 2)) - coords(:, 2) + min(coords(:, 2));
dff = h5read('20180706_Run04.h5', '/Data/dff');
labels = h5read('20180706_Run04.h5', '/Data/labels');
spikes = double(h5read('20180706_Run04.h5', '/Data/spikes'));
bzones = sepRegions(labels, coords, 'plot');

% Convolving spikes:
taud = 0.1;
limexp = 1;
expkern = exp(- (0:0.001:limexp) / taud);
spikesconv = zeros(size(spikes));
for i = 1:size(spikes, 1)
    spikestemp = conv(spikes(i, :), expkern);
    spikesconv(i, :) = spikestemp(1:size(spikes, 2));
end
spikes = spikesconv;
linear = 1;

% Train & test set:
pourc = 0.8;
trainkeep = randperm(size(spikes, 2)-1, round(pourc * size(spikes, 2)));
spikestrain = spikes(:, trainkeep);
yspikestrain = spikes(:, trainkeep+1);
testkeep = find(sum((1:(size(spikes, 2)-1))' == sort(trainkeep), 2) == 0)';
spikestest = spikes(:, testkeep);
yspikestest = spikes(:, testkeep+1);

% Parameters of network:
n0 = size(spikes, 1);
n1 = 100;
W1 = randn(n1, n0) * sqrt(2 / n0);
b1 = zeros(n1, 1);
W2 = randn(n0, n1) * sqrt(2 / n1);
b2 = zeros(n0, 1);
nite = 50;
alpha = 0.05;
J = zeros(1, nite);
Jtest = zeros(1, nite);
modplot = 1;
if exist('linear', 'var') == 0; linear = 0; end
figure
for i = 1:nite
    % Forward for training set:
    Z1 = W1 * spikestrain + b1;
    A1 = Z1 .* (Z1 > 0);
    Z2 = W2 * A1 + b2;
    if linear == 1
        A2 = Z2 .* (Z2 > 0);
        J(i) = mean(sum((yspikestrain - A2).^2));
    else
        A2 = 1 ./ (1 + exp(-Z2));
        J(i) = - mean(sum(yspikestrain .* log(A2) + (1-yspikestrain) .* log(1-A2)));
    end
    % Forward for test set:
    Z1test = W1 * spikestest + b1;
    A1test = Z1test .* (Z1test > 0);
    Z2test = W2 * A1test + b2;
    if linear == 1
        A2test = Z2test .* (Z2test > 0);
        Jtest(i) = mean(sum((yspikestest - A2test).^2));
    else
        A2test = 1 ./ (1 + exp(-Z2test));
        Jtest(i) = - mean(sum(yspikestest .* log(A2test) + (1-yspikestest) .* log(1-A2test)));
    end
    % Plotting:
    if mod(i, modplot) == 0
        pause(0.000001)
        hold off
        plot(J(1:i))
        hold on
        plot(Jtest(1:i))
        legend('Training set cost function', 'Test set cost function')
        xlabel('Iteration', 'Interpreter', 'latex')
    end
    % Backprop:
    if linear == 1
        dA2 = -2*A2 .* (yspikestrain - A2);
        dZ2 = 1 * (dA2 > 0);
    else
        dZ2 = A2 - yspikestrain;
    end
    db2 = sum(dZ2, 2);
    dW2 = dZ2 * A1';
    dZ1 = (W2' * dZ2) .* (1 * (Z1 > 0));
    db1 = sum(dZ1, 2);
    dW1 = dZ1 * spikestrain';
    % Gradient:
    W1 = W1 - alpha*dW1;
    b1 = b1 - alpha*db1;
    W2 = W2 - alpha*dW2;
    b2 = b2 - alpha*db2;
end

% Trying algorithm:
Z1s = W1 * spikes(:, 1:end-1) + b1;
A1s = Z1s .* (Z1s > 0);
Z2s = W2 * A1s + b2;
A2s = 1 ./ (1 + exp(-Z2s));
figure
for i = 1:5
    subplot(5, 1, i)
    hold on
    temp = randperm(size(spikes, 1), 1);
    plot(spikes(temp, :))
    plot([0, A2s(temp, :)])
end
        


% Defining constants:
n = 100;
S = round(rand(n, n))*2 -1;
pattern1 = ones(n, n);
pattern1(floor(0.5*n)+1:n, :) = -1;
pattern2 = ones(n, n);
pattern2(:, floor(0.5*n)+1:n) = -1;
% First plot:
figure
subplot(4, 4, [3, 4, 7, 8])
image(pattern1, 'CDataMapping', 'scaled')
subplot(4, 4, [11, 12, 15, 16])
image(pattern2, 'CDataMapping', 'scaled')
subplot(4, 4, [5, 6, 9, 10])
image(S, 'CDataMapping', 'scaled')
% Hebbian rule for patterns:
W = 0.5 * (pattern1(:)*pattern1(:)' + pattern2(:)*pattern2(:)');
W = W - diag(W);
theta = 0 * ones(n^2, 1);
% Main algorithm:
epoch = 100000;
for i = 1:epoch
    Stemp = W * S(:);
    jrand = randperm(n^2, 1);
    sjtemp = (Stemp(jrand) > theta(jrand))*2 -1;
    S(jrand) = sjtemp;
    if mod(i, 1) == 0
        pause(0.000001)
        subplot(4, 4, [5, 6, 9, 10])
        image(S, 'CDataMapping', 'scaled')
    end
end







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

outij = corClust(valij, 0.99);
maxij = max(outij(2, :));
valijnew = zeros(0, size(valij, 2));



%% Building new signals:

nunew = size(valij, 1) - size(outij, 2) + length(unique(outij(2, :)));
nsignew = valij;
nsignew(outij(1, :), :) = [];
nsignew = [zeros(length(unique(outij(2, :))), size(valij, 2)); nsignew];
for i = 1:length(unique(outij(2, :)))
    outfind = find(outij(2, :) == i);
    nsignew(i, :) = mean(valij(outij(1, outfind), :));
end



%% Doing that for all the zones:

% cormax = 0.65;
% nsignal = cell(size(bzones));
% for i = 1:size(nsignal, 1)
%     for j = 1:size(nsignal, 2)
%         % Information on progress:
%         fprintf('\n\n\n\n\n\n\n\nBrain zone %.0f out of %.0f \n', [size(nsignal, 2)*(i-1)+j, size(nsignal, 1)*size(nsignal, 2)]);
%         % Compute correlations:
%         valij = nvalues(bzones{i, j}, :);
%         outij = corClust(valij, cormax);
%         % New signal values:
%         nsignew = valij;
%         nsignew(outij(1, :), :) = [];
%         nsignew = [zeros(length(unique(outij(2, :))), size(valij, 2)); nsignew];
%         for k = 1:length(unique(outij(2, :)))
%             outfind = find(outij(2, :) == k);
%             nsignew(k, :) = mean(valij(outij(1, outfind), :));
%         end
%         % Allocating to nsignal:
%         nsignal{i, j} = nsignew;
%     end
% end



%% K-means:

addpath('~/Science/Hippolyte/Small-programs')
addpath('~/Science/Hippolyte/Small-programs/braindivision')
   
nite = 30;
rmsestop = 0.001;
nsignal = cell(size(bzones));
for i = 1:size(nsignal, 1)
    for j = 1:size(nsignal, 2)
        % Information on progress:
        fprintf('\n\n\n\n\n\n\n\nBrain zone %.0f out of %.0f \n', [size(nsignal, 2)*(i-1)+j, size(nsignal, 1)*size(nsignal, 2)]);
        % Take values and the number of clusters:
        valij = nvalues(bzones{i, j}, :);
        ncentij = round(size(valij, 1) / 20);
        % New clustering:
        [centij, minclustij, ~] = iKmeans(valij, ncentij, nite, 1, rmsestop);
        clustij = unique(minclustij);
        nsignew = zeros(length(clustij), size(valij, 2));
        for k = 1:length(clustij)
            nsignew(k, :) = mean(valij(minclustij == clustij(k), :));
        end
        % Allocating to nsignal:
        nsignal{i, j} = nsignew;
    end
end

sigtot = [];
for i = 1:size(nsignal, 1)
    for j = 1:size(nsignal, 2)
        sigtot = [sigtot; nsignal{i, j}];
    end
end

% % t = 0:1:(4*pi);
% % s1 = sin(t); s2 = 2*sin(t);
% % sdata = [s1; s2; s2+0.2*randn(size(t)); 10*s1; 0.3*randn(size(t))+0.5; -s1; -s2; -s2+0.05*randn(size(t))];
% % figure; hold on; for i = 1:8; plot(sdata(i, :), '+:'); end
% % [clust, assign, ~] = Kmeans(sdata, 3, nite, 1)
% 
% ss1 = mvnrnd([-6, 6], [1, 0.3; 0.3, 1], 1000);
% ss2 = mvnrnd([3, 3], [1, 0; 0, 2], 1200);
% ss3 = mvnrnd([0, -6], [1, -0.9; -0.9, 1], 950);
% ss4 = [mean(ss1); mean(ss2); mean(ss3)];
% sst = [ss1; ss2; ss3];
% % figure; hold on; plot(sst(:, 1), sst(:, 2), '.')
% [s1, s2, s3, s4] = iKmEM(sst, 3, 3, 0);
% figure
% subplot(1, 2, 1)
% hold on
% [~, s2bis] = max(s2, [], 2);
% for i = 1:size(s1, 1)
%     plot(sst(s2bis == i, 1), sst(s2bis == i, 2), '.', 'Color', rand(1, 3))
% end
% plot(ss4(:, 1), ss4(:, 2), '+g')
% plot(s3(:, 1), s3(:, 2), '+r')
% for i = 1:3
%     plot(s1{i, 2}(1), s1{i, 2}(2), '+k')
% end
% subplot(1, 2, 2)
% plot(s4)
% 
% 
% 
% gg1 = mvnrnd([0, 0], [1, 0.9; 0.9, 1], 1000);
% gg2 = mvnrnd([3, 3], [0.9, 0; 0, 0.9], 600);
% gg = [gg1; gg2];
% figure
% hold on
% plot(gg(:, 1), gg(:, 2), '.')
% [s1, s2, s3] = iKmEM(gg, 2, 0, 1, 0.001);
% % plot(s1(:, 1), s1(:, 2), '+r')
% plot(gg(s2(:, 1)>s2(:, 2), 1), gg(s2(:, 1)>s2(:, 2), 2), '.', 'Color', rand(1, 3))
% plot(gg(s2(:, 1)<=s2(:, 2), 1), gg(s2(:, 1)<=s2(:, 2), 2), '.', 'Color', rand(1, 3))


 









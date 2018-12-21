clear; close all; clc



%% Going to the right address:

sinstim = 0;
switch sinstim
    case 0
        cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-05-24/Run 08/Analysis/HDF5'
    case 1
        cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-06-14/Run 03/Analysis/HDF5'
end




%% We can get the info of the file:

hippo = 'h5hippo.h5';
stim = h5read(hippo, '/Data/Stimulus');
refco = h5read(hippo, '/Data/RefCoordinates');
dff = h5read(hippo, '/Data/Values');
dffn = (dff - mean(dff, 2)) ./ std(dff, [], 2);



%% Based on this info, we can obtain the data inside:

addpath(genpath('~/Science/Hippolyte'))
Params.period = 1;
[out, variables, perinf] = h5stimreg(hippo, Params);
figure
subplot(2, 1, 1)
hold on
[~, n1] = max(out.R2score);
plot(perinf{3}(n1, :))
plot(sum(out.coef(n1, :) .* variables(:, 2:end), 2) + out.intercept(n1))
plot(perinf{2})
title('Highest R2 score signal, with regression and stimulus', 'Interpreter', 'latex')
subplot(2, 1, 2)
hold on
n2 = randperm(size(perinf{3}, 1), 1);
plot(perinf{3}(n2, :))
plot(sum(out.coef(n2, :) .* variables(:, 2:end), 2) + out.intercept(n2))
plot(perinf{2})
title('Random signal, with regression and stimulus', 'Interpreter', 'latex')



%% Plotting neurons with high F-stat and high coefficients:

%Parameters:
nneu = size(dff, 1);
q = 10;
% Coefficients:
coef = sqrt(sum((out.coef.^2), 2));
[maxcoef, indcoef] = sort(coef, 'descend');
indcK = sort(indcoef(1:ceil(nneu/q)));
% F-stat:
fstat = out.Fstat;
[maxfstat, indfstat] = sort(fstat, 'descend');
indfK = sort(indfstat(1:ceil(nneu/q)));
% Taking indices that are in both vectors:
compa = sum(indcK == indfK', 2);
indfin = indcK(compa == 1);
% Computing F-stat and coefficients for neurons kept:
indfinf = ((1 ./ fstat(indfin)) - (1 ./ max(fstat(indfin)))) ./ (1 ./ min(fstat(indfin)));
indfinc = ((1 ./ coef(indfin)) - (1 ./ max(coef(indfin)))) ./ (1 ./ min(coef(indfin)));
indind = (indfinf.*indfinc - min(indfinf.*indfinc)) ./ max(indfinf.*indfinc);

% Plotting:
% load('MaskDatabase.mat')
% brain = [275, 1, 94, 114, 260];
% width = 621; height = 1406; Zs = 138;
% x = 1:width;
% y = height:-1:1;
% z = 1:Zs;
% [X, Y, Z] = meshgrid(x, y, z);
figure
hold on
grid on
axis equal
% for i = 1:length(brain)
%     masktemp = full(MaskDatabaseOutlines(:, brain(i)));
%     sparsekeep = find(masktemp == 1);
%     scatter3(X(sparsekeep), Y(sparsekeep), Z(sparsekeep), 10, 'filled', 'MarkerFaceAlpha', 1/16, 'MarkerFaceColor', [0.8, 0.8, 0.8])
% end
scatter3(refco(indfin, 1), refco(indfin, 2), refco(indfin, 3), [], [indind, indind, indind], '.')
axis equal
title('Neurons with high F-stat and high regression coefficients', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')

figure
subplot(1, 2, 1)
hold on
grid on
axis equal
scatter3(refco(indfin, 1), refco(indfin, 2), refco(indfin, 3), [], [indfinf, 0.9*ones(size(indfinf)), 0.9*ones(size(indfinf))], '.')
axis equal
title('Neurons with high F-stat', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')
subplot(1, 2, 2)
hold on
grid on
axis equal
scatter3(refco(indfin, 1), refco(indfin, 2), refco(indfin, 3), [], [0.9*ones(size(indfinf)), indfinc, 0.9*ones(size(indfinf))], '.')
axis equal
title('Neurons with high regression coefficients', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')



%% DBscanning:

refcok = refco(indfin, :);
Params = struct;
Params.dmax = 0.015;
Params.minpt = 20;
% Params.weight = 1 ./ (indind + 1);
% Params.pweight = 10;
Params.plot = 0;
Params.fprint = 500;
Params.monitor = 0;
out = dbscan(refcok, Params);

indref1 = find(refcok(:, 1) <= 0.25);
indref2 = find(refcok(:, 1) > 0.25);
refco1 = refcok(indref1, :);
refco2 = refcok(indref2, :);
clust1 = dbscan(refco1, Params);
clust2 = dbscan(refco2, Params);
clust2(clust2 ~= -1) = clust2(clust2 ~= -1) + max(clust1);
clust = [clust1; clust2];
figure
hold on
axis equal
view(0, 90)
grid on
title('DBSCAN clustering on 3D data', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
clustmax = max(clust);
for i = 1:clustmax
    if i <= max(clust1)
        scatter3(refco1(clust1 == i, 1), refco1(clust1 == i, 2), refco1(clust1 == i, 3), ...
                 50, 'filled', 'MarkerFaceAlpha', 5/8, 'MarkerFaceColor', rand(1, 3))
    else
        scatter3(refco2(clust2 == i, 1), refco2(clust2 == i, 2), refco2(clust2 == i, 3), ...
                 50, 'filled', 'MarkerFaceAlpha', 5/8, 'MarkerFaceColor', rand(1, 3))
    end        
end



%% Testing h5regress function:

[out2, variables2, perinf2] = h5regress(hippo, variables(:, 2:end)');
isequal(out, out2) % nope, because signal renormalized













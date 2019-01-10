clear; close all; clc



%% Going to the right address:

sinstim = 0;
switch sinstim
    case 0
        cd '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/2018-05-24/Run 08/Analysis/HDF5'
    case 1
        cd '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2018-06-14/Run 03/Analysis/HDF5'
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
outout = dbscan(refcok, Params);

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

% [out2, variables2, perinf2] = h5regress(hippo, variables(:, 2:end)');
% isequal(out, out2) % nope, because signal renormalized



% datatest = rand(1000, 2);
% figure; plot(datatest(:, 1), datatest(:, 2), '.'); grid on
% Params = struct;
% Params.dmax = 0.05;
% Params.minpt = 10;
% weight = rand(1000, 1);
% weight(weight > 0.9) = 100;
% Params.weight = weight;
% Params.pweight = 1;
% Params.plot = 1;
% Params.fprint = 500;
% Params.monitor = 0;
% out = dbscan(datatest, Params);



%% Plotting depending on coefficients:

coeff = out.coef;
nstim = size(coeff, 2);
keepneu = cell(1, nstim);
figure
for i = 1:nstim
    cotemp = coeff(:, i);
    keepneu{i} = find(cotemp < mean(cotemp)-2*std(cotemp) | cotemp > mean(cotemp)+2*std(cotemp));
    keepcoef = coeff(keepneu{i}, i);
    keepcoef = (keepcoef - min(keepcoef)) ./ (max(keepcoef) - min(keepcoef));
    subplot(2, 2, i)
    scatter3(refco(keepneu{i}, 1), refco(keepneu{i}, 2), refco(keepneu{i}, 3), ...
                 [], [keepcoef, 1-keepcoef, keepcoef], '.')
    axis equal
    view(0, 90)
    grid on
    title('Coefficient to one of 4 stimuli', 'Interpreter', 'latex')
    xlabel('x-axis', 'Interpreter', 'latex')
    ylabel('y-axis', 'Interpreter', 'latex')
    zlabel('z-axis', 'Interpreter', 'latex')
end


% dffnb = dffn(bzones{1, 1}, :);
% meand = zeros(size(dffnb, 1));
% for i = 2:size(dffnb, 1)
%     for j = 1:i
%         meand(i, j) = sqrt(sum((dffnb(i, :)-dffnb(j, :)).^2));
%     end
% end
% meand = meand + meand';
% figure
% image(meand, 'CDataMapping', 'scaled')
% colorbar
        


%% Comparing vestibular and thermotaxis:

Params.period = 0;
[outV, variablesV, perinfV] = h5stimreg(hippo, Params);
addpath('/home/ljp/Science/Guillaume/Thermotaxis/Datasets')
% hippoT = '20171115_Run06Tset=27.h5';
hippoT = '20180111_Run04Tset=14.h5';
[outT, variablesT, perinfT] = h5stimreg(hippoT, Params);
figure
subplot(2, 1, 1)
hold on
[~, n1] = max(outT.R2score);
plot(perinfT{3}(n1, :))
plot(sum(outT.coef(n1, :) .* variablesT(:, 2:end), 2) + outT.intercept(n1))
plot(perinfT{2})
title('Highest R2 score signal, with regression and stimulus', 'Interpreter', 'latex')
subplot(2, 1, 2)
hold on
n2 = randperm(size(perinfT{3}, 1), 1);
plot(perinfT{3}(n2, :))
plot(sum(outT.coef(n2, :) .* variablesT(:, 2:end), 2) + outT.intercept(n2))
plot(perinfT{2})
title('Random signal, with regression and stimulus', 'Interpreter', 'latex')

% Plotting neurons with high F-stat and high coefficients:
%Parameters:
dffT = h5read(hippoT, '/Data/dff');
nneuT = size(dffT, 1);
qT = 10;
% Coefficients:
coefT = sqrt(sum((outT.coef.^2), 2));
[maxcoefT, indcoefT] = sort(coefT, 'descend');
indcKT = sort(indcoefT(1:ceil(nneuT/qT)));
% F-stat:
fstatT = outT.Fstat;
[maxfstatT, indfstatT] = sort(fstatT, 'descend');
indfKT = sort(indfstatT(1:ceil(nneuT/qT)));
% Taking indices that are in both vectors:
compaT = sum(indcKT == indfKT', 2);
indfinT = indcKT(compaT == 1);
% Computing F-stat and coefficients for neurons kept:
indfinfT = ((1 ./ fstatT(indfinT)) - (1 ./ max(fstatT(indfinT)))) ./ (1 ./ min(fstatT(indfinT)));
indfincT = ((1 ./ coefT(indfinT)) - (1 ./ max(coefT(indfinT)))) ./ (1 ./ min(coefT(indfinT)));
indindT = (indfinfT.*indfincT - min(indfinfT.*indfincT)) ./ max(indfinfT.*indfincT);

% Plotting:
refcoT = double(h5read(hippoT, '/Data/zbrain_coords'));
refcoT(:, 2) = max(refcoT(:, 2)) - refcoT(:, 2);
refcoT = ((refcoT - min(refcoT)) ./ (max(refcoT) - min(refcoT))) .* (max(refco) - min(refco)) + min(refco);
figure; subplot(1, 2, 1); scatter3(refco(:, 1), refco(:, 2), refco(:, 3)); axis equal; subplot(1, 2, 2); scatter3(refcoT(:, 1), refcoT(:, 2), refcoT(:, 3)); axis equal
figure
hold on
grid on
axis equal
scatter3(refcoT(indfinT, 1), refcoT(indfinT, 2), refcoT(indfinT, 3), [], [indindT, indindT, indindT], '.')
axis equal
title('Neurons with high F-stat and high regression coefficients (thermotaxis)', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')

figure
subplot(1, 2, 1)
hold on
grid on
axis equal
scatter3(refcoT(indfinT, 1), refcoT(indfinT, 2), refcoT(indfinT, 3), [], [indfinfT, 0.9*ones(size(indfinfT)), 0.9*ones(size(indfinfT))], '.')
axis equal
title('Neurons with high F-stat (thermotaxis)', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')
subplot(1, 2, 2)
hold on
grid on
axis equal
scatter3(refcoT(indfinT, 1), refcoT(indfinT, 2), refcoT(indfinT, 3), [], [0.9*ones(size(indfinfT)), indfincT, 0.9*ones(size(indfinfT))], '.')
axis equal
title('Neurons with high regression coefficients (thermotaxis)', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')















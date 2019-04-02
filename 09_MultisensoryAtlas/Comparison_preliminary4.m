% clear; close all; clc
% addpath(genpath('/home/ljp/Programs'))
% addpath(genpath('/home/ljp/Science/Hippolyte'))
close all; clc
% 
% 
% 
% %% Loading datasets:
% 
% load('zgridGeo01.mat')
% load('zgridGui01.mat')



%% Assiging new datasets:

zgrid = zgridGeo + zgridGui;
% zgridGeoStep = choseSubset(zgridGeo, 'step');
% zgridGeoSine = choseSubset(zgridGeo, 'sine');
% zgridGuiHot = choseSubset(zgridGui, 'hot');
% zgridGuiCold = choseSubset(zgridGui, 'cold');
% zgridGuiNeut = choseSubset(zgridGui, 'neutral');



%% Studying correlation distibutions:

corre = cell(5, 1);
segments = ["step", "sine", "hot", "cold", "neutral"];
for i = 1:5
    zgridtemp = choseSubset(zgrid, segments(i));
    lzgrid = length(zgridtemp);
    corretemp = zeros(lzgrid, 2);
    for j = 1:lzgrid
        corretemp(j, 1) = mean(zgridtemp.Zcorvect{j});
        corretemp(j, 2) = var(zgridtemp.Zcorvect{j});
    end
    corre{i} = corretemp;
end

cormat = zeros(5, 4);
for i = 1:5
    cormat(i, 1) = mean(corre{i}(:, 1));
    cormat(i, 2) = var(corre{i}(:, 1)); 
    cormat(i, 3) = mean(corre{i}(:, 2));
    cormat(i, 4) = var(corre{i}(:, 2)); 
end
titles = ["Mean correlation across datasets", "Mean correlation variance across datasets", ...
          "Correlation variance across datasets", "Correlation variance variance across datasets"];
subsets = ["Vestibular step", "Vestibular sine", "Thermotaxis hot", "Thermotaxis cold", "Thermotaxis neutral"];
figure
for i = 1:4
    subplot(4, 1, i)
    cortemp = cormat(:, i);
    plot(cortemp, ':o', 'MarkerFaceColor', [0, 0, 0])
    title(titles(i), 'Interpreter', 'latex')
    difftemp = max(cortemp) - min(cortemp);
    axis([0.9, 5.1, min(cortemp)-0.1*difftemp, max(cortemp)+0.1*difftemp])
    xticks(1:5)
    xticklabels(subsets)
    grid on
end



%% Flattening:

zgridGeoStepf = flatten(choseSubset(zgridGeo, 'step'));
zgridGeoSinef = flatten(choseSubset(zgridGeo, 'sine'));
zgridGuiHotf = flatten(choseSubset(zgridGui, 'hot'));
zgridGuiColdf = flatten(choseSubset(zgridGui, 'cold'));
zgridGuiNeutf = flatten(choseSubset(zgridGui, 'neutral'));
    

temp = zgridGuiHot(1);
cor = temp.Zcorvect{1};
dff = h5read(char(temp.paths), '/Data/Brain/Analysis/DFF');
stim = h5read(char(temp.paths), '/Data/Stimulus/RandomPulses/Trace');
stim = reshape(stim, 1, length(stim));
dffn = (dff - mean(dff, 2)) ./ std(dff, [], 2);
stimn = (stim - mean(stim)) / std(stim);
cor2 = sum(dffn .* stimn, 2) ./ sqrt(sum(dffn.^2, 2) .* sum(stimn.^2, 2));
figure
for i = 1:5
    subplot(6, 1, i)
    ntemp = randperm(size(dff, 1), 1);
    plot(dff(ntemp, :))
    title(['Neuron ', num2str(ntemp), ', with correlation ', num2str(cor(ntemp)), ', and new correlation ', num2str(cor2(ntemp))], 'Interpreter', 'latex')
end
subplot(6, 1, 6)
plot(stim)



%% Problem with the correlations: have to do it again...

zgrid_copy = zgrid;
lz = length(zgrid);
for i = 1:lz
    % Defining object:
    ztemp = zgrid(i);
    % Getting dff:
    dff = h5read(char(ztemp.paths), '/Data/Brain/Analysis/DFF');
    dffn = (dff - mean(dff, 2)) ./ std(dff, [], 2);
    % Getting stimulus:
    findthermo = strfind(ztemp.comments, 'thermotaxis');
    if isempty(findthermo)
        stim = h5read(char(ztemp.paths), '/Data/Stimulus/vestibular1/motorAngle');
    else
        stim = h5read(char(ztemp.paths), '/Data/Stimulus/RandomPulses/Trace');
    end
    stim = reshape(stim, 1, length(stim));
    stimn = (stim - mean(stim)) / std(stim);
    % Computing correlation:
    cor = sum(dffn .* stimn, 2) ./ sqrt(sum(dffn.^2, 2) .* sum(stimn.^2, 2));
    % Modifying zgrid object:
    zgrid.Zcorvect{i} = cor;
    [x, y, z, d] = size(zgrid.Zcorrelations);
    Zcortemp = zeros(x, y, z);
    indtemp = find(zgrid.Zcorrelations(:, :, :, i) ~= 0);
    for t = 1:length(indtemp)
        Ttemp = zgrid.Zneurons(:, :, :, i);
        Ttemp = Ttemp(indtemp(t));
        Zcortemp(t) = mean(cor(Ttemp{1}));
    end
    zgrid.Zcorrelations(:, :, :, i) = Zcortemp;
    fprintf('Iteration %.0f out of %.0f. \n', [i, lz]);
end



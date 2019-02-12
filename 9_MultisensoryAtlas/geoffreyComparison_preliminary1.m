clear; close all; clc
addpath(genpath('/home/ljp/Programs'));
addpath(genpath('/home/ljp/Science/Hippolyte/small-programs'))



%% Loading 3 step-stimulus files from Geoffrey:

h1 = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run08';
h2 = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run11';
h3 = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run15';
H = {h1, h2, h3};
nh = length(H);
names = cell(size(H));
for i = 1:nh
    names{i} = H{i}(end-14:end);
end



%% Getting data from files:

% Stimulus:
stim = cell(size(H));
for i = 1:nh
    stim{i} = h5read(H{i}, '/Data/Stimulus/vestibular1/motorAngle');
%     % Problem with 200 first points for run 11, 2018-05-24:
%     if string(H{i}) == '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run11'
%         stim{i} = stim{i}(201:end);
%     end
end
% DFF:
dff = cell(size(H));
for i = 1:nh
    dff{i} = h5read(H{i}, '/Data/Brain/Analysis/DFF');
%     % Problem with 200 first points for run 11, 2018-05-24:
%     if string(H{i}) == '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run11'
%         dff{i} = dff{i}(:, 201:end);
%     end
end
% Coordinates:
coord = cell(size(H));
for i = 1:nh
    coord{i} = h5read(H{i}, '/Data/Brain/ZBrainCoordinates');
end



%% Normalizing data:

% Stimulus:
stimn = cell(size(H));
for i = 1:nh
    stimn{i} = (stim{i} - mean(stim{i})) / std(stim{i});
end
% DFF:
dffn = cell(size(H));
for i = 1:nh
    dffn{i} = (dff{i} - mean(dff{i}, 2)) ./ std(dff{i}, [], 2);
end
% Plotting random neurons dff:
figure
for i = 1:5
    subplot(5, 1, i)
    hold on
    for j = 1:nh
        plot(dffn{j}(randperm(size(dffn{j}, 1), 1), :))
    end
end



%% Computing correlation coefficients:

method = "Correlation analysis";
cor = cell(size(H));
for i = 1:nh
    cor{i} = sum(dffn{i} .* stimn{i}, 2) ./ sqrt(sum(dffn{i}.^2, 2) .* sum(stimn{i}.^2, 2));
end



%% Plotting all correlation coefficients to have an idea:

figure
for i = 1:nh
    subplot(nh, 4, ((i-1)*4+1):i*4-1)
    plot(cor{i}, '.')
    title(['Correlation coefficients for dataset ', num2str(i)], 'Interpreter', 'latex')
    xlabel('Neuron', 'Interpreter', 'latex')
end
subplot(nh, 4, 4:4:(4*nh))
vartemp = zeros(nh, 1);
for i = 1:nh
    cortemp = cor{i};
    cortemp(isnan(cortemp)) = [];
    vartemp(i) = var(cortemp);
end
image(vartemp, 'CDataMapping', 'scaled'); colorbar
title('Correlation coefficient variance for all datasets', 'Interpreter', 'latex')

% %% Loading reference brain:
% 
% [X, meta] = nrrdread('/home/ljp/Science/GeoffreysComputer/RefBrains/zBrain_Elavl3-H2BRFP_198layers.nhdr');



%% Creating a linspace map for the brain and assigning neurons:

% Taking extrema:
extmax = zeros(nh, 3);
extmin = zeros(nh, 3);
for i = 1:nh
    for j = 1:3
        extmax(i, j) = max(coord{i}(:, j));
        extmin(i, j) = min(coord{i}(:, j));
    end
end
extmax = max(extmax);
extmin = min(extmin);
extrema = [extmin; extmax];

% % Defining increment for grid:
% increment = 0.005;
% 
% % Creating grid:
% xgrid = extrema(1, 1):increment:extrema(2, 1); 
% lx = length(xgrid);
% ygrid = extrema(1, 2):increment:extrema(2, 2);
% ly = length(ygrid);
% zgrid = extrema(1, 3):increment:extrema(2, 3);
% lz = length(zgrid);
% Tgrid = cell(lx-1, ly-1, lz-1, nh);
% Cgrid = zeros(lx-1, ly-1, lz-1, nh);
% % Filling grid with neurons:
% tic
% for ix = 1:(lx-1)
%     for iy = 1:(ly-1)
%         for iz = 1:(lz-1)
%             for id = 1:nh
%                 % For each dataset:
%                 xtemp = (xgrid(ix) <= coord{id}(:, 1) & coord{id}(:, 1) < xgrid(ix+1));
%                 ytemp = (ygrid(iy) <= coord{id}(:, 2) & coord{id}(:, 2) < ygrid(iy+1));
%                 ztemp = (zgrid(iz) <= coord{id}(:, 3) & coord{id}(:, 3) < zgrid(iz+1));
%                 Ttemp = find(xtemp & ytemp & ztemp);
%                 Tgrid{ix, iy, iz, id} = Ttemp;
%                 Cgrid(ix, iy, iz, id) = mean(cor{id}(Ttemp));
%             end
%         end
%     end
%     fprintf('%.2f %% done, in %.3f seconds. \n', [100*ix/(lx-1), toc]);
% end
% % Building structure:
% gridStruct = struct;
% gridStruct.method = method;
% gridStruct.names = names;
% gridStruct.xgrid = xgrid;
% gridStruct.ygrid = ygrid;
% gridStruct.zgrid = zgrid;
% gridStruct.Tgrid = Tgrid;
% gridStruct.Cgrid = Cgrid;

% Loading file for increment 0.005:
load('/home/ljp/Science/Hippolyte/gridStruct.mat')
method = gridStruct.method;
names = gridStruct.names;
xgrid = gridStruct.xgrid;
ygrid = gridStruct.ygrid;
zgrid = gridStruct.zgrid;
Tgrid = gridStruct.Tgrid;
Cgrid = gridStruct.Cgrid;
increment = xgrid(2) - xgrid(1);

% Computing mean correlation:
Cgridm = mean(Cgrid, 4);



%% Plotting:

% Building meshgrid:
xmesh = (xgrid(1:end-1)+xgrid(2:end)) / 2;
ymesh = (ygrid(1:end-1)+ygrid(2:end)) / 2;
zmesh = (zgrid(1:end-1)+zgrid(2:end)) / 2;
[Xmesh, Ymesh, Zmesh] = meshgrid(xmesh, ymesh, zmesh);
Xmesh = permute(Xmesh, [2, 1, 3]);
Ymesh = permute(Ymesh, [2, 1, 3]);
Zmesh = permute(Zmesh, [2, 1, 3]);

% Grid coordinates:
g_coord = [Xmesh(:), Ymesh(:), Zmesh(:)];

% Getting rid of unnecessary neurons:
ridvalues = [-0.1, 0.1];
[Cgridmlayed, g_coord] = ridNeurons(Cgridm(:), ridvalues, g_coord);

% Plotting:
Ccolor = corr2col(Cgridmlayed);
figure
scatter3(g_coord(:, 1), g_coord(:, 2), g_coord(:, 3), [], Ccolor, 'filled')
axis equal
title('Grid inside Zbrain with correlation', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')



%% Plotting histogram of correlations:

figure
hist(Cgridm(:), 50)
title('Histogram of correlations averaged on grid', 'Interpreter', 'latex')
xlabel('Correlation', 'Interpreter', 'latex')
ylabel('Number of redundancies', 'Interpreter', 'latex')
xticks((-0.5:0.1:0.5))



%% Comparing original zbrain with grid:

% Grid coordinates:
g_coord = [Xmesh(:), Ymesh(:), Zmesh(:)];

% Interval:
ridvalues = [-0.15, 0.15];

% Plotting:
figure 
for i = 1:3
    subplot(3, 2, (i-1)*2+1)
    [Ctemp, coordtemp] = ridNeurons(cor{i}, ridvalues, coord{i});
    scatter3(coordtemp(:, 1), coordtemp(:, 2), coordtemp(:, 3), [], corr2col(Ctemp), 'filled')
    axis equal
    title(['Original neurons and correlations from dataset ', num2str(i)], 'Interpreter', 'latex')
    xlabel('x-axis', 'Interpreter', 'latex')
    ylabel('y-axis', 'Interpreter', 'latex')
    zlabel('z-axis', 'Interpreter', 'latex')
    subplot(3, 2, 2*i)
    Cgridtemp = Cgrid(:, :, :, i);
    [Ctemp, coordtemp] = ridNeurons(Cgridtemp(:), ridvalues, g_coord);
    scatter3(coordtemp(:, 1), coordtemp(:, 2), coordtemp(:, 3), [], corr2col(Ctemp), 'filled')
    axis equal
    title(['Grid neurons and correlations from dataset ', num2str(i)], 'Interpreter', 'latex')
    xlabel('x-axis', 'Interpreter', 'latex')
    ylabel('y-axis', 'Interpreter', 'latex')
    zlabel('z-axis', 'Interpreter', 'latex')
end



%% Trying new functions create_gridStruct and add_to_gridStruct:

method = 'Correlation analysis';
increment = 0.005;
for i = 1:length(H)
    % Loading data:
    stimtemp = h5read(H{i}, '/Data/Stimulus/vestibular1/motorAngle');
    dfftemp = h5read(H{i}, '/Data/Brain/Analysis/DFF');
    coordtemp = h5read(H{i}, '/Data/Brain/ZBrainCoordinates');
    % Normalizing:
    stimntemp = (stimtemp - mean(stimtemp)) / std(stimtemp);
    dffntemp = (dfftemp - mean(dfftemp, 2)) ./ std(dfftemp, [], 2);
    % Correlation coefficient:
    cortemp = sum(dffntemp .* stimntemp, 2) ./ sqrt(sum(dffntemp.^2, 2) .* sum(stimntemp.^2, 2));
    % Computing gridStruct:
    if i == 1
        % If first example:
        gridStruct = create_gridStruct(method, names{i}, coordtemp, cortemp, increment);
    else
        % Otherwise adding info:
        gridStruct = add_to_gridStruct(gridStruct, names{i}, coordtemp, cortemp);
    end
end
% Plotting with meshgrid:
xgrid = gridStruct.xgrid;
ygrid = gridStruct.ygrid;
zgrid = gridStruct.zgrid;
Tgrid = gridStruct.Tgrid;
Cgrid = gridStruct.Cgrid;
Cgridm = mean(Cgrid, 4);
% Building meshgrid:
xmesh = (xgrid(1:end-1)+xgrid(2:end)) / 2;
ymesh = (ygrid(1:end-1)+ygrid(2:end)) / 2;
zmesh = (zgrid(1:end-1)+zgrid(2:end)) / 2;
[Xmesh, Ymesh, Zmesh] = meshgrid(xmesh, ymesh, zmesh);
Xmesh = permute(Xmesh, [2, 1, 3]);
Ymesh = permute(Ymesh, [2, 1, 3]);
Zmesh = permute(Zmesh, [2, 1, 3]);
% Grid coordinates:
g_coord = [Xmesh(:), Ymesh(:), Zmesh(:)];
% Getting rid of unnecessary neurons:
ridvalues = [-0.15, 0.15];
[Cgridmlayed, g_coord] = ridNeurons(Cgridm(:), ridvalues, g_coord);
% Plotting:
Ccolor = corr2col(Cgridmlayed);
figure
scatter3(g_coord(:, 1), g_coord(:, 2), g_coord(:, 3), [], Ccolor, 'filled')
axis equal
title('Grid inside Zbrain with correlation', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')



%% Comparing the two gridStruct:

% gridStruct:
load('/home/ljp/Science/Hippolyte/gridStruct.mat')
method = gridStruct.method;
names = gridStruct.names;
xgrid = gridStruct.xgrid;
ygrid = gridStruct.ygrid;
zgrid = gridStruct.zgrid;
Tgrid = gridStruct.Tgrid;
Cgrid = gridStruct.Cgrid;
increment = xgrid(2) - xgrid(1);
% gridStruct2:
load('/home/ljp/Science/Hippolyte/gridStruct2.mat')
method2 = gridStruct2.method;
names2 = gridStruct2.names;
xgrid2 = gridStruct2.xgrid;
ygrid2 = gridStruct2.ygrid;
zgrid2 = gridStruct2.zgrid;
Tgrid2 = gridStruct2.Tgrid;
Cgrid2 = gridStruct2.Cgrid;
increment2 = xgrid2(2) - xgrid2(1);



%% Checking what is wrong:

n1 = randperm(length(xgrid)-1, 1);
n2 = randperm(length(ygrid)-1, 1);
n3 = randperm(length(zgrid)-1, 1);
X = [xgrid(n1), xgrid(n1+1)], Y = [ygrid(n2), ygrid(n2+1)], Z = [zgrid(n3), zgrid(n3+1)]
set = 1;
{coord{set}(Tgrid{n1, n2, n3, set}, :), cor{set}(Tgrid{n1, n2, n3, set}, :), Cgrid(n1, n2, n3, set)}
set = 2;
{coord{set}(Tgrid{n1, n2, n3, set}, :), cor{set}(Tgrid{n1, n2, n3, set}, :), Cgrid(n1, n2, n3, set)}
set = 3;
{coord{set}(Tgrid{n1, n2, n3, set}, :), cor{set}(Tgrid{n1, n2, n3, set}, :), Cgrid(n1, n2, n3, set)}
% coord{set}(Tgrid2{n1, n2, n3, set}, :)
% cor{set}(Tgrid2{n1, n2, n3, set}, :)
% Cgrid2(n1, n2, n3, set)









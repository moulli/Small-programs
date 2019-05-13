clear; close all; clc
load('/home/ljp/Science/Hippolyte/MaskDatabase.mat')



%% Test:

figure
hold on
pixize = [0.000798, 0.000798, 0.002];
for i = [1, 78, 94, 114, 260, 275]
    temp = full(MaskDatabase(:, i));
    temp = reshape(temp, height, width, Zs);
    temp = find(temp == 1);
    [yt, xt, zt] = ind2sub([height, width, Zs], temp);
    xt = pixize(1) * xt - (pixize(1) / 2);
    yt = pixize(2) * yt - (pixize(2) / 2);
    zt = pixize(3) * zt - (pixize(3) / 2);
    scatter3(xt, yt, zt, '.')
end
axis equal
grid on


%% Length of each brain part:

lenpart = zeros(length(MaskDatabaseNames), 1);
for i = 1:length(MaskDatabaseNames)
    temp = full(MaskDatabase(:, i));
    temp = reshape(temp, height, width, Zs);
    temp = find(temp == 1);
    lenpart(i) = length(temp);
end
figure
plot(lenpart, '.')
grid on


%% Building zlabels005:

% Defining zgrid:
method = 'ZBraingrid object stacking labels associated to a 0.005mm grid object';
zbrainsize = [0.496, 1.122, 0.276];
increment = 0.005;
gridsize = floor(zbrainsize ./ increment);
orientation = 'RAS';
zlabels005 = ZBraingrid(method, gridsize, orientation);

% Adding labels:
pixize = [0.000798, 0.000798, 0.002];
for i = 1:length(MaskDatabaseNames)
    % Defining structure:
    stemp = struct;
    stemp.name = MaskDatabaseNames{i};
    stemp.path = '/home/ljp/Science/Hippolyte/MaskDatabase.mat';
    stemp.orientation = 'RPS';
    % Computing coordinates:
    temp = full(MaskDatabase(:, i));
    temp = reshape(temp, height, width, Zs);
    temp = find(temp == 1);
    [yt, xt, zt] = ind2sub([height, width, Zs], temp);
    xt = pixize(1) * xt - (pixize(1) / 2);
    yt = pixize(2) * yt - (pixize(2) / 2);
    zt = pixize(3) * zt - (pixize(3) / 2);
    stemp.coordinates = [xt, yt, zt];
    stemp.correlation = ones(length(xt), 1);
    % Add to ZBraingrid object:
    addDataset(zlabels005, stemp);
    disp(zlabels005)
end
    
    




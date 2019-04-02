clear; close all; clc


%% This file to learn how to manipulate HDF5 files on Matlab


%% Going to the right address:

cd '~/Science/GeoffreysComputer/RLS/Data/2018-05-24/Run 08/Analysis/HDF5'
% cd '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2018-05-24/Run 08/Analysis/HDF5'


%% We can get the info of the file:

hippo = 'h5hippo.h5';
h5disp(hippo)
h5info(hippo)


%% Based on this info, we can obtain the data inside:

coordinates = h5read(hippo, '/Data/RefCoordinates');
% size(coordinates)
values = h5read(hippo, '/Data/Values');
% size(values)
stimulus = h5read(hippo, '/Data/Stimulus');
% size(stimulus)
df = h5read(hippo, '/Data/df_aligned')';
% size(df)


%% Normalizing the data:

nstimulus = (stimulus-mean(stimulus)) ./ std(stimulus);
nvalues = (values-mean(values, 2)) ./ std(values, [], 2);


%% Computing covariance between neurons and stimulus:

covariance = sum(nvalues.*nstimulus, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimulus.^2, 2));
figure
hist(covariance, 50)
grid on
title('Histogram of covariances between activations and stimulus (centered on 0: symetry)', 'Interpreter', 'latex')
xlabel('Covariance', 'Interpreter', 'latex')
ylabel('Number of apparitions', 'Interpreter', 'latex')
xticks([-1, -0.5, 0, 0.5, 1])


%% Plotting a nice representation:

figure
scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), [], [0.2, 0.2, 0.2], '.')
hold on
% scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 150, 'filled', 'MarkerFaceAlpha', 3/16, 'MarkerFaceColor', [0.8, 0.8, 0.8])
axis equal
title('All recorded neurons in the brain', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')

HnegC = find(covariance <= -0.65);
Hncoord = coordinates(HnegC, :);
MnegC = find(-0.65 < covariance & covariance <= -0.5);
Mncoord = coordinates(MnegC, :);
LnegC = find(-0.5 < covariance & covariance <= -0.45);
Lncoord = coordinates(LnegC, :);
HposC = find(covariance >= 0.65);
Hpcoord = coordinates(HposC, :);
MposC = find(0.5 <= covariance & covariance < 0.65);
Mpcoord = coordinates(MposC, :);
LposC = find(0.45 <= covariance & covariance < 0.5);
Lpcoord = coordinates(LposC, :);
noC = find(-0.45 < covariance & covariance < 0.45);
nocoord = coordinates(noC, :);

figure
scatter3(Hncoord(:, 1), Hncoord(:, 2), Hncoord(:, 3), [], [1, 0, 0], '.')
hold on
scatter3(Mncoord(:, 1), Mncoord(:, 2), Mncoord(:, 3), [], 	[1, 0.7, 0.3], '.')
scatter3(Lncoord(:, 1), Lncoord(:, 2), Lncoord(:, 3), [], 	[1, 0.9, 0.5], '.')
scatter3(Hpcoord(:, 1), Hpcoord(:, 2), Hpcoord(:, 3), [], [0, 0, 1], '.')
scatter3(Mpcoord(:, 1), Mpcoord(:, 2), Mpcoord(:, 3), [], [0.3, 0.7, 1], '.')
scatter3(Lpcoord(:, 1), Lpcoord(:, 2), Lpcoord(:, 3), [], [0.5, 0.9, 1], '.')
hold off
legend('High negative correlation', 'Medium negative correlation', 'Low negative correlation', 'High positive correlation', ...
       'Medium positive correlation', 'Low positive correlation')
axis equal
title('Correlation between neurons activation and stimulus', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')


%% What if we only keep the positive/negative part of the stimulus:

stimpos = nstimulus .* (nstimulus >= 0);
covpos = sum(nvalues.*stimpos, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(stimpos.^2, 2));
stimneg = nstimulus .* (nstimulus <= 0);
covneg = sum(nvalues.*stimneg, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(stimneg.^2, 2));

negCpos = find(covpos <= -0.5);
npcoord = coordinates(negCpos, :);
posCpos = find(covpos >= 0.5);
ppcoord = coordinates(posCpos, :);
negCneg = find(covneg <= -0.5);
nncoord = coordinates(negCneg, :);
posCneg = find(covneg >= 0.5);
pncoord = coordinates(posCneg, :);

figure
scatter3(npcoord(:, 1), npcoord(:, 2), npcoord(:, 3), [], [0.4, 0, 0.4], 'o')
hold on
scatter3(ppcoord(:, 1), ppcoord(:, 2), ppcoord(:, 3), [], 	[1, 0, 1], '.')
scatter3(nncoord(:, 1), nncoord(:, 2), nncoord(:, 3), [], [0, 0.4, 0], '.')
scatter3(pncoord(:, 1), pncoord(:, 2), pncoord(:, 3), [], [0, 1, 0], 'o')
hold off
legend('Negative correlation to positive stimulus', 'Positive correlation to positive stimulus', 'Negative correlation to negative stimulus (phase shifting)', ...
       'Positive correlation to negative stimulus')
axis equal
title('High and medium correlation between neurons activation and positive or negative stimulus', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')


%% Plotting with these new indications:

HnegC = find(covneg <= -0.65);
Hncoord = coordinates(HnegC, :);
MnegC = find(-0.65 < covneg & covneg <= -0.5);
Mncoord = coordinates(MnegC, :);
LnegC = find(-0.5 < covneg & covneg <= -0.45);
Lncoord = coordinates(LnegC, :);
lnegC = find(-0.45 < covneg & covneg <= -0.4);
lncoord = coordinates(lnegC, :);
HposC = find(covpos >= 0.65);
Hpcoord = coordinates(HposC, :);
MposC = find(0.5 <= covpos & covpos < 0.65);
Mpcoord = coordinates(MposC, :);
LposC = find(0.45 <= covpos & covpos < 0.5);
Lpcoord = coordinates(LposC, :);
lposC = find(0.4 <= covpos & covpos < 0.45);
lpcoord = coordinates(lposC, :);

figure
scatter3(Hncoord(:, 1), Hncoord(:, 2), Hncoord(:, 3), [], [1, 0, 0], '.')
hold on
scatter3(Mncoord(:, 1), Mncoord(:, 2), Mncoord(:, 3), [], 	[1, 0.7, 0.3], '.')
scatter3(Lncoord(:, 1), Lncoord(:, 2), Lncoord(:, 3), [], 	[1, 0.9, 0.5], '.')
scatter3(lncoord(:, 1), lncoord(:, 2), lncoord(:, 3), [], 	[1, 1, 0.6], '.')
scatter3(Hpcoord(:, 1), Hpcoord(:, 2), Hpcoord(:, 3), [], [0, 0, 1], '.')
scatter3(Mpcoord(:, 1), Mpcoord(:, 2), Mpcoord(:, 3), [], [0.3, 0.7, 1], '.')
scatter3(Lpcoord(:, 1), Lpcoord(:, 2), Lpcoord(:, 3), [], [0.5, 0.9, 1], '.')
scatter3(lpcoord(:, 1), lpcoord(:, 2), lpcoord(:, 3), [], [0.7, 1, 1], '.')
hold off
legend('Very high correlation to negative stimulus', 'High correlation to negative stimulus', 'Medium correlation to negative stimulus', 'Low correlation to negative stimulus', ...
       'Very high correlation to positive stimulus', 'High correlation to positive stimulus', 'Medium correlation to positive stimulus', 'Low correlation to positive stimulus')
axis equal
title('Correlation between neurons activation and positive or negative stimulus', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')



%% Which neurons correlated to the stimulus are correlated to each other:

% keep = 53;
% [maxi, indi] = sort(covariance);
% indix1 = indi(1:keep);
% indix2 = indi(end-keep+1:end);
indix1 = posCpos;
indix2 = negCneg;
indix = [indix1; indix2];
objetx = nvalues(indix, :);
covx = cov(objetx');
figure
image(covx, 'CDataMapping', 'scaled')
colorbar
title('Covariance between neurons for positive correlation to positive signal, and negative correlation to negative signal', 'Interpreter', 'latex')


% Keeping neurons that activate:
figure
subplot(1, 2, 1)
covx1 = cov(nvalues(indix1, :)');
mcovx1 = mean(covx1);
[~, indd1] = sort(mcovx1);
len1 = length(indix1);
lin1 = [linspace(0, 1, len1)', linspace(0, 1, len1)', linspace(0, 0, len1)'];
scatter3(coordinates(indix1(indd1), 1), coordinates(indix1(indd1), 2), coordinates(indix1(indd1), 3), [], lin1, '.')
axis equal
title('Correlation between neurons of posive correlation to positive stimulus (black low, yellow high)', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')
subplot(1, 2, 2)
covx2 = cov(nvalues(indix2, :)');
mcovx2 = mean(covx2);
[~, indd2] = sort(mcovx2);
len2 = length(indix2);
lin2 = [linspace(0, 1, len2)', linspace(0, 1, len2)', linspace(0, 0, len2)'];
scatter3(coordinates(indix2(indd2), 1), coordinates(indix2(indd2), 2), coordinates(indix2(indd2), 3), [], lin2, '.')
axis equal
title('Correlation between neurons of negative correlation to negative stimulus (black low, yellow high)', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')

% thrld = 0.85;
% actstim1 = find(mean(covx1) > thrld);
% actstim2 = find(mean(covx2) > thrld);
% actstimn1 = find(mean(covx1) <= thrld);
% actstimn2 = find(mean(covx2) <= thrld);
% figure
% scatter3(coordinates(indix1(actstim1), 1), coordinates(indix1(actstim1), 2), coordinates(indix1(actstim1), 3), [], [1, 1, 0], '.')
% axis equal
% hold on
% scatter3(coordinates(indix2(actstim2), 1), coordinates(indix2(actstim2), 2), coordinates(indix2(actstim2), 3), [], [1, 1, 0], '.')
% scatter3(coordinates(indix1(actstimn1), 1), coordinates(indix1(actstimn1), 2), coordinates(indix1(actstimn1), 3), [], [0.6, 0.6, 0.6], '.')
% scatter3(coordinates(indix2(actstimn2), 1), coordinates(indix2(actstimn2), 2), coordinates(indix2(actstimn2), 3), [], [0.6, 0.6, 0.6], '.')


%% Finding the clusters:

thrld = 0.95;

lenx1 = length(covx1);
cupx1 = triu(covx1, 1);
cupx1bin = find(cupx1 > thrld);
X1bin = mod(cupx1bin, lenx1) + lenx1*(mod(cupx1bin, lenx1)==0);
Y1bin = ceil(cupx1bin/lenx1);
lenx1bin = length(cupx1bin);
clusttemp = [];
for k = 1:lenx1bin
    if isempty(find(clusttemp == X1bin(k), 1)) == 0 && isempty(find(clusttemp == Y1bin(k), 1)) == 0
        colX1 = ceil(find(clusttemp == X1bin(k))/size(clusttemp, 1));
        colY1 = ceil(find(clusttemp == Y1bin(k))/size(clusttemp, 1));
        if colX1 == colY1
            continue
        else
            indtemp = (clusttemp(:, colY1) ~= 0);
            clusttemp(indtemp, colX1) = clusttemp(indtemp, colY1);
            clusttemp = clusttemp(:, [1:colY1-1, colY1+1:end]);
        end
    elseif isempty(find(clusttemp == X1bin(k), 1)) == 0
        colx1 = ceil(find(clusttemp == X1bin(k))/size(clusttemp, 1));
        clusttemp = [clusttemp; zeros(1, size(clusttemp, 2))];
        clusttemp(end, colx1) = Y1bin(k);
    elseif isempty(find(clusttemp == Y1bin(k), 1)) == 0
        colx1 = ceil(find(clusttemp == Y1bin(k))/size(clusttemp, 1));
        clusttemp = [clusttemp; zeros(1, size(clusttemp, 2))];
        clusttemp(end, colx1) = X1bin(k);
    else
        clusttemp = [clusttemp, zeros(size(clusttemp, 1), 1)];
        clusttemp = [clusttemp; zeros(2, size(clusttemp, 2))];
        clusttemp(end-1:end, end) = [X1bin(k); Y1bin(k)];
    end
end
clust1 = [];
for k = 1:size(clusttemp, 2)
    if nnz(clusttemp(:, k)) > 2
        clust1 = [clust1, clusttemp(:, k)];
    end
end


%% Periodicity:

stimwind = nstimulus(49:102);
lenwind = length(stimwind);
lenstim = length(nstimulus);
duree = lenstim - lenwind + 1;
covstim = zeros(duree, 1);
for i = 1:duree
    covstim(i) = sum(stimwind.*nstimulus(i:i+lenwind-1), 2) ./ sqrt(sum(stimwind.^2, 2).*sum(nstimulus(i:i+lenwind-1).^2, 2));
end
[~, covmax] = sort(covstim);
covmax = covmax(end-14:end);
covmax = sort(covmax);
period = covmax(2) - covmax(1);
minper = period / 2;

lencovperiod = length(stimulus) / period;
covperiod = zeros(length(indix), lencovperiod);
for i = 1:lencovperiod
    covperiod(:, i) = sum(nvalues(indix, 1:period).*nvalues(indix, period*(i-1)+1:i*period), 2) ./ sqrt(sum(nvalues(indix, 1:period).^2, 2).*sum(nvalues(indix, period*(i-1)+1:i*period).^2, 2));
end
covperiod = mean(covperiod, 2);

lencovminper = length(stimulus) / minper;
covminper = zeros(length(indix), lencovminper);
for i = 1:lencovminper
    covminper(:, i) = sum(nvalues(indix, 1:minper).*nvalues(indix, minper*(i-1)+1:i*minper), 2) ./ sqrt(sum(nvalues(indix, 1:minper).^2, 2).*sum(nvalues(indix, minper*(i-1)+1:i*minper).^2, 2));
end
covminper = mean(covminper, 2);

figure
plot(covperiod, '.')
hold on
plot(covminper, '.')
legend('Period of positive and negative stimulus', 'Period of only positive or negative stimulus')
title('Periodicity for the neurons the most correlated to stimulus', 'Interpreter', 'latex')
xlabel('Neuron', 'Interpreter', 'latex')
ylabel('Periodicity of signal', 'Interpreter', 'latex')
grid on



Hper = find(covperiod >= 0.65);
Hpcoord = coordinates(indix(Hper), :);
Mper = find(0.65 > covperiod & covperiod >= 0.5);
Mpcoord = coordinates(indix(Mper), :);
Lper = find(0.5 > covperiod & covperiod >= 0.35);
Lpcoord = coordinates(indix(Lper), :);
lper = find(0.35 > covperiod);
lpcoord = coordinates(indix(lper), :);

figure
scatter3(Hpcoord(:, 1), Hpcoord(:, 2), Hpcoord(:, 3), [], [0, 0, 0], '.')
hold on
scatter3(Mpcoord(:, 1), Mpcoord(:, 2), Mpcoord(:, 3), [], 	[0.3, 0.3, 0.3], '.')
scatter3(Lpcoord(:, 1), Lpcoord(:, 2), Lpcoord(:, 3), [], 	[0.6, 0.6, 0.6], '.')
scatter3(lpcoord(:, 1), lpcoord(:, 2), lpcoord(:, 3), [], [0.8, 0.8, 0.8], '.')
hold off
legend('Very high periodicity', 'High periodicity', 'Medium periodicity', 'Low periodicity')
axis equal
title('Periodicity of neurons signals', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')

        
        
%% Clustering:

data = coordinates(indix, :);
% Parameters:
maxdist = 0.025;
minpt = 20;
coreffect = 10;
plotclust = 1;
% data = [1+0.1*randn(100, 1), 1+0.1*randn(100, 1), 1+0.1*randn(100, 1);
%               -1+0.1*randn(100, 1), 0.1*randn(100, 1), -2+0.1*randn(100, 1);
%               randn(30, 1), 0.5+0.5*randn(30, 1), -0.5+1.5*randn(30, 1)];
clusters = dbscan3(data, maxdist, minpt, covariance(indix), coreffect, plotclust);
clustmax = max(clusters);
sizeclust = zeros(clustmax, 1);
for i = 1:clustmax
    sizeclust(i) = sum(clusters == i);
end
figure
plot(sizeclust, 'x-.')
hold on
axis([1, clustmax, 0, max(sizeclust)])
title('Number of neurons per cluster', 'Interpreter', 'latex')
xlabel('Cluster', 'Interpreter', 'latex')
ylabel('Number of neurons', 'Interpreter', 'latex')
grid on


%% Clustering, separating positive and negative correlation:

% Computing clusters:
data1 = coordinates(indix1, :);
maxdist1 = 0.02;
minpt1 = 35;
coreffect1 = 1000;
clusters1 = dbscan3(data1, maxdist1, minpt1, covariance(indix1), coreffect1);
data2 = coordinates(indix2, :);
maxdist2 = 0.02;
minpt2 = 35;
coreffect2 = 1000;
clusters2 = dbscan3(data2, maxdist2, minpt2, covariance(indix2), coreffect2);
clusters2 = clusters2 + max(clusters1);
clusters2(clusters2 == max(clusters1)-1) = -1;
clusters = [clusters1; clusters2];
clustmax = max(clusters);

% Colors we will use:
colmin = 0.1;
colmax = 0.9;
numcol = ceil(clustmax/6)+1;
coloruse = [linspace(colmax, colmax, numcol)', linspace(colmin, colmax, numcol)', linspace(colmin, colmin, numcol)';
            linspace(colmax, colmin, numcol)', linspace(colmax, colmax, numcol)', linspace(colmin, colmin, numcol)';
            linspace(colmin, colmin, numcol)', linspace(colmax, colmax, numcol)', linspace(colmin, colmax, numcol)';
            linspace(colmin, colmin, numcol)', linspace(colmax, colmin, numcol)', linspace(colmax, colmax, numcol)';
            linspace(colmin, colmax, numcol)', linspace(colmin, colmin, numcol)', linspace(colmax, colmax, numcol)';
            linspace(colmax, colmax, numcol)', linspace(colmin, colmin, numcol)', linspace(colmax, colmin, numcol)'];
getridcol = numcol:numcol:6*numcol;
coloruse(getridcol, :) = [];
permcolor = randperm(size(coloruse, 1));
        
% Plotting the new clusters:
figure
hold on
axis equal
view(0, 90)
grid on
title('DBSCAN clustering on 3D data, separating positive and negative correlation', 'Interpreter', 'latex')
xlabel('x-axis', 'Interpreter', 'latex')
ylabel('y-axis', 'Interpreter', 'latex')
zlabel('z-axis', 'Interpreter', 'latex')
data = coordinates(indix, :);
for i = 1:clustmax
    scatter3(data(clusters == i, 1), data(clusters == i, 2), data(clusters == i, 3), ...
             50, 'filled', 'MarkerFaceAlpha', 5/8, 'MarkerFaceColor', coloruse(permcolor(i), :))
end









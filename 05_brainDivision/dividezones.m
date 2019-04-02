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
   
% nite = 30;
% rmsestop = 0.001;
% nsignal = cell(size(bzones));
% for i = 1:size(nsignal, 1)
%     for j = 1:size(nsignal, 2)
%         % Information on progress:
%         fprintf('\n\n\n\n\n\n\n\nBrain zone %.0f out of %.0f \n', [size(nsignal, 2)*(i-1)+j, size(nsignal, 1)*size(nsignal, 2)]);
%         % Take values and the number of clusters:
%         valij = nvalues(bzones{i, j}, :);
%         ncentij = round(size(valij, 1) / 20);
%         % New clustering:
%         [centij, minclustij, ~] = iKmeans(valij, ncentij, nite, 1, rmsestop);
%         clustij = unique(minclustij);
%         nsignew = zeros(length(clustij), size(valij, 2));
%         for k = 1:length(clustij)
%             nsignew(k, :) = mean(valij(minclustij == clustij(k), :));
%         end
%         % Allocating to nsignal:
%         nsignal{i, j} = nsignew;
%     end
% end
addpath('~/Science/Hippolyte')
load('nsigiKd.mat')

% nite = 30;
% rmsestop = 0.001;
% nsigiKd = cell(size(bzones));
% for i = 1:size(nsigiKd, 1)
%     for j = 1:size(nsigiKd, 2)
%         % Information on progress:
%         fprintf('\n\n\n\n\n\n\n\nBrain zone %.0f out of %.0f \n', [size(nsigiKd, 2)*(i-1)+j, size(nsigiKd, 1)*size(nsigiKd, 2)]);
%         % Take values and the number of clusters:
%         valij = nvalues(bzones{i, j}, :);
%         ncentij = round(size(valij, 1) / 20);
%         % New clustering:
%         [a, b, c, d] = iKmeans(valij, ncentij, nite, 1, rmsestop);
%         nsigiKd{i, j} = {a, b, c, d};
%     end
% end
load('nsigiKd2.mat')
nsigiKd2f = cell(size(bzones));
for i = 1:size(nsigiKd2f, 1)
    for j = 1:size(nsigiKd2f, 2)
        centij = nsigiKd2{i, j}{1};
        minclustij = nsigiKd2{i, j}{2};
        clustij = unique(minclustij);
        nsignew = zeros(length(clustij), size(valij, 2));
        valij = nvalues(bzones{i, j}, :);
        for k = 1:length(clustij)
            nsignew(k, :) = mean(valij(minclustij == clustij(k), :));
        end
        % Allocating to nsignal:
        nsigiKd2f{i, j} = nsignew;
    end
end


% Reorganizing based on brain areas, from anterior to posterior
sortnsig = [5, 1, 2, 3, 4]';
nsigiKd2f = nsigiKd2f(sortnsig, :);

sigtot = [];
for i = 1:size(nsigiKd2f, 1)
    for j = 1:size(nsigiKd2f, 2)
        sigtot = [sigtot; nsigiKd2f{i, j}];
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



%% Deconvolution using BSD:

addpath(genpath('~/Programs/BSD'))

Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = size(sigtot, 2); % Number of time frames.
f = 2.5; % Frequency of acquisition.
Oalg.dt = 1 / f; % interval duration. (s)
Oalg.nNeurons = 1; % Number of neurons.
Oalg.adaptive = 1; % Adaptive. Will use the provided values in P as initializer, estimate the unprovided ones, then iteratively refine the values.
Oalg.iterations = 10; % Maximal number of iterations. Default: 5.
t = Oalg.dt * (0:(length(sigtot(1, :))-1));

tic
[Nspikes, Conv, Palg, Pphys, Oalg] = BSD(sigtot(1, :), Oalg);
toc
figure
subplot(4, 1, 1)
hold on
plot(t, sigtot(1, :), 'Color', [0.5, 0.5, 0.5])
plot(t, Conv, 'green')
legend('Calcium signal', 'Deconvolved signal')
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Fluorescence', 'Interpreter', 'latex')
subplot(4, 1, 2)
plot(t, Nspikes)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Number of spikes in interval', 'Interpreter', 'latex')
subplot(4, 1, 3)
Nbins = Nspikes > Pphys.threshold;
plot(t, Nbins, 'red')
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Number of spikes in interval, approximated', 'Interpreter', 'latex')
% Checking if Nspikes = 2 means that there is two spikes:
expt = exp(- t ./ Pphys.tauDecay) - exp(- t ./ Pphys.tauRise);
trac = zeros(Oalg.Time);
for i = 1:Oalg.Time
    trac(i:end, i) = expt(1:end-i+1)';
end
C=filter(1,[Palg.eta -Palg.eta*(Palg.gamma+Palg.delta),Palg.eta*Palg.delta],Nbins) + Palg.b;
C2=filter(1/Palg.eta,[1, -(Palg.gamma+Palg.delta), Palg.delta],Nbins) + Palg.b;
% Plotting:
subplot(4, 1, 4)
hold on
plot(t, sigtot(1, :), 'Color', [0.5, 0.5, 0.5])
plot(t, C, 'blue')
legend('Calcium signal', 'Approximated signal with whole spikes')
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Fluorescence', 'Interpreter', 'latex')


%% Trying to come up with a way to transform into real spikes:

Nnew = zeros(1, 2*Oalg.Time);
Nstemp = Nspikes;
nmin = 0.5;
nmax = 1.5;
for i = 1:Oalg.Time
    switch i
        case Oalg.Time
            if nmin <= Nstemp(i) && Nstemp(i) <= nmax
                Nnew(2*i-1) = 1;
            elseif nmax < Nstemp(i)
                Nnew(2*i-1) = 1;
                Nnew(2*i) = 1;
            end
        otherwise
            if nmin <= Nstemp(i) && Nstemp(i) <= nmax
                Nnew(2*i-1) = 1;
            elseif nmax < Nstemp(i)
                Nnew(2*i-1) = 1;
                Nstemp(2*i) = 1;
            elseif Nstemp(i) < nmin && nmin <= (Nstemp(i) + Nstemp(i+1)) && (Nstemp(i) + Nstemp(i+1)) <= nmax
                Nnew(2*i) = 1;
                Nstemp(i+1) = 0;
            end
    end
end
Cnew=filter(1,[Palg.eta -Palg.eta*(Palg.gamma+Palg.delta),Palg.eta*Palg.delta],Nnew) + Palg.b;

gg1 = Nnew;
gg2 = reshape([Nspikes; zeros(1, size(Nspikes, 2))], [1, length(gg1)]);
C1=filter(1,[Palg.eta -Palg.eta*(Palg.gamma+Palg.delta),Palg.eta*Palg.delta],gg1) + Palg.b;
C2=filter(1,[Palg.eta -Palg.eta*(Palg.gamma+Palg.delta),Palg.eta*Palg.delta],gg2) + Palg.b;
figure
subplot(2, 1, 1)
hold on
plot(0:0.2:(0.4*(Oalg.Time-1)+0.2), gg1)
plot(0:0.2:(0.4*(Oalg.Time-1)+0.2), gg2)
subplot(2, 1, 2)
hold on
plot(0:0.2:(0.4*(Oalg.Time-1)+0.2), C1)
plot(0:0.2:(0.4*(Oalg.Time-1)+0.2), C2)

figure
hold on
plot(t, sigtot(1, :), 'Color', [0.5, 0.5, 0.5])
plot(0:0.2:(0.4*(Oalg.Time-1)+0.2), Cnew, 'green')
plot(0:0.2:(0.4*(Oalg.Time-1)+0.2), Nnew, 'red')
legend('Calcium signal', 'Deconvolved signal approximation', 'Spikes approximation')
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Fluorescence', 'Interpreter', 'latex')        
% Really long, let's go with values here.


%% Superesolution:

% res = 10;
% Oalg.superResolution = res;
% tic
% [Nspikessup, Convsup, Palgsup, Pphyssup, Oalgsup] = BSD(sigtot(1, :), Oalg);
% toc
% figure
% subplot(2, 1, 1)
% hold on
% plot(1:Oalg.Time, Nspikes)
% plot(1/res:1/res:Oalg.Time, Nspikessup)
% subplot(2, 1, 2)
% hold on
% Nspikesres = reshape([Nspikes; zeros(res-1, size(Nspikes, 2))], [1, length(Nspikessup)]);
% plot(1/res:1/res:Oalg.Time, filter(1,[Palg.eta -Palg.eta*(Palg.gamma+Palg.delta),Palg.eta*Palg.delta], Nspikesres) + Palg.b)
% plot(1/res:1/res:Oalg.Time, filter(1,[Palg.eta -Palg.eta*(Palg.gamma+Palg.delta),Palg.eta*Palg.delta], Nspikessup) + Palg.b)
% Nbinsup = Nspikessup > Pphyssup.threshold;
% plot(1/res:1/res:Oalg.Time, filter(1,[Palg.eta -Palg.eta*(Palg.gamma+Palg.delta),Palg.eta*Palg.delta], Nbinsup) + Palg.b)





%% PARAMETERS:
% SPIKES OR NUMBER OF SPIKES PER TIME STEP:
spikes = 0;
deconv = 1;
keepone = 0;



%% Deconvolving everything:

Oalg = struct; % Struct of experimental conditions & decoding options.
Oalg.Time = size(sigtot, 2); % Number of time frames.
f = 2.5; % Frequency of acquisition.
Oalg.dt = 1 / f; % interval duration. (s)
Oalg.nNeurons = size(sigtot, 1); % Number of neurons.
Oalg.adaptive = 1; % Adaptive. Will use the provided values in P as initializer, estimate the unprovided ones, then iteratively refine the values.
Oalg.iterations = 10; % Maximal number of iterations. Default: 5.
t = Oalg.dt * (0:(length(sigtot(1, :))-1));

% [Nspikes, Conv, Palg, Pphys, Oalg] = BSD(sigtot, Oalg);
% Deconv = {Nspikes, Conv, Palg, Pphys, Oalg};
load('Deconv.mat')
[Nspikes, Conv, Palg, Pphys, Oalg] = Deconv{:};
if deconv == 1
    Nspikes = Conv - Palg.b'; % USING DECONVOLVED SIGNAL INSTEAD OF SPIKES
end


%% Adding stimulus to data:


% if spikes == 1
%     Nspikes = double(Nspikes > Pphys.threshold');
% end
% 
% stimn = (stimulus - mean(stimulus)) ./ std(stimulus);
% stimv = [0, stimulus(2:end)-stimulus(1:end-1)];
% stimvn = (stimv - mean(stimv)) ./ std(stimv);
% stima = [0, stimv(2:end)-stimv(1:end-1)];
% stiman = (stima - mean(stima)) ./ std(stima);
% if spikes == 1
%     Nspikesn = Nspikes;
% else
%     Nspikesn = (Nspikes) ./ std(Nspikes(:));
% end
% 
% if keepone == 1
%     Nstot = [ones(1, Oalg.Time); stimn; stimvn; stiman; Nspikesn];
% else
%     Nstot = [stimn; stimvn; stiman; Nspikesn];
% end
% nones = size(Nstot, 1) - size(Nspikes, 1);

if spikes == 1
    Nspikes = double(Nspikes > Pphys.threshold');
end

dt = 0.4;
expkern =  exp(-(0:dt:(20*2.16))/2.16);

stimn = (stimulus - mean(stimulus)) ./ std(stimulus);
stimnp = stimn .* (stimn > 0);
stimnn = stimn .* (stimn < 0);
ksp = convInd(stimnp, expkern, 1, length(stimulus));
ksn = convInd(stimnn, expkern, 1, length(stimulus));

stimv = gradient(stimn);
stimvn = (stimv - mean(stimv)) ./ std(stimv);
stimvnp = stimvn .* (stimvn > 0);
stimvnn = stimvn .* (stimvn < 0);
kvp = convInd(stimvnp, expkern, 1, length(stimulus));
kvn = convInd(stimvnn, expkern, 1, length(stimulus));

stima = gradient(stimvn);
stiman = (stima - mean(stima)) ./ std(stima);
stimanp = stiman .* (stiman > 0);
stimann = stiman .* (stiman < 0);
kap = convInd(stimanp, expkern, 1, length(stimulus));
kan = convInd(stimann, expkern, 1, length(stimulus));
if spikes == 1
    Nspikesn = Nspikes;
else
    Nspikesn = (Nspikes) ./ std(Nspikes(:));
end

if keepone == 1
    Nstot = [ones(1, Oalg.Time); ksp; ksn; kvp; kvn; kap; kan; Nspikesn];
else
    Nstot = [ksp; ksn; kvp; kvn; kap; kan; Nspikesn];
end
nones = size(Nstot, 1) - size(Nspikes, 1);

%% Finding weights:

% Training and test set:
trainperc = 0.85;
Ntrain = Nstot(:, 1:ceil(trainperc*Oalg.Time));
Ntest = Nstot(:, ceil(trainperc*Oalg.Time)+1:end);


% Initialization:
Wopt = 0.012 * (rand(Oalg.nNeurons, size(Nstot, 1)) - 0.5);
% Optimization algorithm:
lambda = 10;
mup = 1.2;
mum = 0.5;
dmax = 50;
Delta0 = 0.01 * ones(size(Wopt));
Loss{1} = zeros(size(Wopt));
Loss{2} = zeros(size(Wopt));
dgrad = 0.0005;
nite = 500;
rmsei = zeros(nite, 1);
figure
for i = 1:nite
    % Initialization : 
    datain = Ntrain(:, 1:end-1);
    [dinx, diny] = size(datain);
    dataout = Ntrain(nones+1:end, 2:end);
    % Compute output:
    ytemp = Wopt * datain;
    if spikes == 1
        yi = 1 ./ (1 + exp(-ytemp));
    else
        yi = ytemp .* (ytemp > 0);
    end
    % Compute RMSE:
    rmsei(i) = sqrt(mean(mean((yi - dataout).^2)));
    % Compute gradient:
    errori = yi - dataout;
    if keepone == 1
        Wreg = [zeros(size(Wopt, 1), 1), Wopt(:, 2:end)];
    else
        Wreg = Wopt;
    end       
    gradJ = (1 / diny) * (errori * datain' + lambda * Wreg);
%     % Resilient backprop update:
%     Loss{1} = Loss{2};
%     Loss{2} = gradJ;
%     losstemp = Loss{1} .* Loss{2};
%     Delta = mup*Delta0.*(losstemp > 0) + mum*Delta0.*(losstemp < 0) ...
%                         + Delta0.*(losstemp == 0);
%     Delta(Delta > dmax) = dmax;
%     Delta0 = Delta;
%     DeltaW = -Delta0.*(Loss{2} > 0) + Delta0.*(Loss{2} < 0);
    % Gradient descent update:
    DeltaW = - dgrad * gradJ;
    % Weights update:
    Wopt = Wopt + DeltaW;
    % Info on progress:
    if mod(i, 10) == 0
        fprintf('Iteration %.0f out of %.0f \n', [i, nite])
    end
    % Plotting rmse:
    pause(0.000001)
    plot(1:i, rmsei(1:i))
end
% figure
% plot(rmsei)
% figure; image(Loss{2}, 'CDataMapping', 'scaled'); colorbar
% figure; image(errori, 'CDataMapping', 'scaled'); colorbar
% figure; image(datain', 'CDataMapping', 'scaled'); colorbar
% figure; image(DeltaW, 'CDataMapping', 'scaled'); colorbar
% figure; hold on; plot(yi(1, :)- dataout(1, :))
figure; hold on; plot(yi(1, :)); plot(dataout(1, :))
% figure; plot(sum(errori, 2))
% figure; plot(sort(sum(errori, 2)))
    
% Test set:
datatest = Ntest(:, 1:end-1);
[dtex, dtey] = size(datatest);
dataoutest = Ntest(nones+1:end, 2:end);
% Compute output:
ytemp = Wopt * datatest;
yit = ytemp .* (ytemp > 0);
% Compute RMSE:
rmsetest = sqrt(mean(mean((yit - dataoutest).^2)));
figure; plot(mean(dataoutest), '.'); hold on; plot(mean(yit), '.'); plot(mean((yit - dataoutest).^2), '.')
% figure; hold on; plot(yi(1, :)- dataoutest(1, :))
figure; hold on; plot(yit(1, :)); plot(dataoutest(1, :))
fprintf('RMSE for training set is %.4f, for test set is %.4f \n', [rmsei(end), rmsetest]);

if keepone == 1
    col = ['k'; 'b'; 'b'; 'r'; 'r'; 'y'; 'y'];
else
    col = ['b'; 'b'; 'r'; 'r'; 'y'; 'y'];
end
figure
for i = 1:nones
    subplot(nones, 1, i)
    plot(Wopt(:, i), col(i))
    axis([1, size(Wopt, 1), -0.03, 0.03])
end

tickslol = [0.001, 0.002, 0.005];
ticks3d = permute(tickslol, [1, 3, 2]);
Woptbispos = Wopt > ticks3d;
Woptbisneg = Wopt < -ticks3d;
Woptbis = sum(Woptbispos, 3) - sum(Woptbisneg, 3);
xline = [0.5, 0.5];
xinit = 0.5;
xend = size(Woptbis);
figure
hold on
image(Woptbis, 'CDataMapping', 'scaled')
tickslab = num2str([-fliplr(tickslol), 0, tickslol]');
colorbar('Ticks', (-length(tickslol):length(tickslol)),...
         'TickLabels', mat2cell(tickslab, ones(1, size(tickslab, 1)), size(tickslab, 2)))
xline(1) = xline(1) + nones;
line([xline(1), xline(1)], [xinit, xend(1)], 'Color', 'black')
xzoned = zeros(size(nsigiKd2f, 1) * size(nsigiKd2f, 2), 1);
for i = 1:size(nsigiKd2f, 1)
    for j = 1:size(nsigiKd2f, 2)
        if i == size(nsigiKd2f, 1) && j == size(nsigiKd2f, 2)
            xzoned(size(nsigiKd2f, 2)*(i-1)+j) = size(nsigiKd2f{i, j}, 1);
            continue
        end
        xzoned(size(nsigiKd2f, 2)*(i-1)+j) = size(nsigiKd2f{i, j}, 1);
        xline(1) = xline(1) + size(nsigiKd2f{i, j}, 1);
        xline(2) = xline(2) + size(nsigiKd2f{i, j}, 1);
        line([xline(1), xline(1)], [xinit, xend(1)], 'Color', 'black')
        line([xinit, xend(2)], [xline(2), xline(2)], 'Color', 'black')
    end
end
axis([xinit, xend(2), xinit, xend(1)])



%% Mean and std:

xzones1 = cumsum([0; xzoned]);
xzones2 = cumsum([0; nones; xzoned]);
xz1max = length(xzones1) - 1;
xz2max = length(xzones2) - 1;
Woptmean = zeros(xz1max, xz2max);
Woptstd = zeros(xz1max, xz2max);
Woptmax = zeros(xz1max, xz2max);
Woptmin = zeros(xz1max, xz2max);
for i = 1:xz1max
    for j = 1:xz2max
        mattemp = Wopt(xzones1(i)+1:xzones1(i+1), xzones2(j)+1:xzones2(j+1));
        Woptmean(i, j) = mean(mattemp(:));
        Woptstd(i, j) = std(mattemp(:));
        Woptmax(i, j) = max(mattemp(:));
        Woptmin(i, j) = min(mattemp(:));
    end
end
figure
subplot(2, 2, 1)
image(Woptmean, 'CDataMapping', 'scaled')
colorbar
title('Mean synaptic weight for all major areas', 'Interpreter', 'latex')
subplot(2, 2, 2)
image(Woptstd, 'CDataMapping', 'scaled')
colorbar
title('Std of synaptic weight for all major areas', 'Interpreter', 'latex')
subplot(2, 2, 3)
image(Woptmax, 'CDataMapping', 'scaled')
colorbar
title('Max synaptic weight for all major areas', 'Interpreter', 'latex')
subplot(2, 2, 4)
image(Woptmin, 'CDataMapping', 'scaled')
colorbar
title('Min synaptic weight for all major areas', 'Interpreter', 'latex')



%% Mean and std without moelle épinière:

xzones1 = xzones1(1:end-2);
xzones2 = xzones2(1:end-2);
xz1max = length(xzones1) - 1;
xz2max = length(xzones2) - 1;
Woptmean = zeros(xz1max, xz2max);
Woptstd = zeros(xz1max, xz2max);
Woptmax = zeros(xz1max, xz2max);
Woptmin = zeros(xz1max, xz2max);
for i = 1:xz1max
    for j = 1:xz2max
        mattemp = Wopt(xzones1(i)+1:xzones1(i+1), xzones2(j)+1:xzones2(j+1));
        Woptmean(i, j) = mean(mattemp(:));
        Woptstd(i, j) = std(mattemp(:));
        Woptmax(i, j) = max(mattemp(:));
        Woptmin(i, j) = min(mattemp(:));
    end
end
figure
subplot(2, 2, 1)
image(Woptmean, 'CDataMapping', 'scaled')
colorbar
title('Mean synaptic weight for all major brain areas', 'Interpreter', 'latex')
subplot(2, 2, 2)
image(Woptstd, 'CDataMapping', 'scaled')
colorbar
title('Std of synaptic weight for all major brain areas', 'Interpreter', 'latex')
subplot(2, 2, 3)
image(Woptmax, 'CDataMapping', 'scaled')
colorbar
title('Max synaptic weight for all major brain areas', 'Interpreter', 'latex')
subplot(2, 2, 4)
image(Woptmin, 'CDataMapping', 'scaled')
colorbar
title('Min synaptic weight for all major brain areas', 'Interpreter', 'latex')



%% Plotting in continuous:

matstim = [ksp; ksn; kvp; kvn; kap; kan];
ntime = length(ksp);
nclust = size(Nspikesn, 1);
randneu = randperm(nclust, 1);
NspikesnB = zeros(size(Nspikesn));
NspikesnB(:, 1) = Nspikesn(:, 1);
figure
subplot(2, 1, 1)
plot(0, stimulus(1), 'Color', [0.8, 0.8, 0])
subplot(2, 1, 2)
plot(0, Nspikesn(randneu, 1))
hold on
plot(0, NspikesnB(randneu, 1), 'r')
hold off
for i = 2:ntime
    pause(0.001)
    Nstemp = [matstim(:, i-1); NspikesnB(:, i-1)];
    Nstemp = Wopt * Nstemp;
    NspikesnB(:, i) = Nstemp .* (Nstemp > 0);
    subplot(2, 1, 1)
    plot(0:0.4:(0.4*(i-1)), stimulus(1:i), 'Color', [0.8, 0.8, 0])
    subplot(2, 1, 2)
    plot(0:0.4:(0.4*(i-1)), Nspikesn(randneu, 1:i))
    hold on
    plot(0:0.4:(0.4*(i-1)), NspikesnB(randneu, 1:i), 'r')
    hold off
end
subplot(2, 1, 1)
plot(0:0.4:(0.4*(i-1)), stimulus(1:i), 'Color', [0.8, 0.8, 0])
title('Stimulus against time', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')
subplot(2, 1, 2)
plot(0:0.4:(0.4*(i-1)), Nspikesn(randneu, 1:i))
hold on
plot(0:0.4:(0.4*(i-1)), NspikesnB(randneu, 1:i), 'r')
hold off
legend('True signal', 'Approximated signal')
title('Signals against time', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

diff = sqrt(mean((Nspikesn - NspikesnB).^2));
figure
plot(diff, '.')



 








function outkern = convInd(u, v, i1, i2)

%% Function that performs convolution and keep only indicated indices.
%
%
%% Parameters:
%  --u: first vector to convolve.
%  --v: second vector to convolve.
%  --i1: index to which we start to keep information.
%  --i2: index to which we stop to keep information.
%
%
%% Output:
%  --outkern: convolved vector, with i2-i1+1 components.



    %% Initialization:
    
    if ~isrow(u) && ~iscolumn(u) 
        error('Please provide first input as a vector')
    elseif ~isrow(v) && ~iscolumn(v)
        error('Please provide second input as a vector')
    elseif i1 < 1 || i1 > i2 || i2 > length(u)+length(v)-1
        error('Inputs must be as follow: 0 < input3 <= input4 < length(input1)+length(input2)')
    end
    
    
    %% Main code:
    
    outkern = conv(u, v);
    outkern = outkern(i1:i2);   



end
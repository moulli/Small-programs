clear; close all; clc



%% Going to the right address:

sine = 1;
if sine == 1
    cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-06-14/Run 03/Analysis/HDF5'
else
    cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-05-24/Run 08/Analysis/HDF5'
end



%% We can get the info of the file:

hippo = 'h5hippo.h5';



%% Based on this info, we can obtain the data inside:

coordinates = h5read(hippo, '/Data/RefCoordinates');
values = h5read(hippo, '/Data/Values');
stimulus = h5read(hippo, '/Data/Stimulus');
labels = h5read(hippo, '/Data/Labels');



%% Normalizing the data:

nstimulus = (stimulus-mean(stimulus)) ./ std(stimulus);
nvalues = (values-mean(values, 2)) ./ std(values, [], 2);



%% Computing covariance between neurons and stimulus:

covariance = sum(nvalues.*nstimulus, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimulus.^2, 2));



%% Separating positive and negative stimuli:

nstimabs = abs(nstimulus);
nstimpos = nstimulus .* (nstimulus > 0);
nstimneg = -nstimulus .* (nstimulus < 0);

covabs = sum(nvalues.*nstimabs, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimabs.^2, 2));
covpos = sum(nvalues.*nstimpos, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimpos.^2, 2));
covneg = sum(nvalues.*nstimneg, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nstimneg.^2, 2));



%% Normalizing data:

covp = (covpos - mean(covpos)) ./ std(covpos);
covn = (covneg - mean(covneg)) ./ std(covneg);



%% Plotting most correlated neurons:

infp = 1;
supp = 400;
[sortp, indp] = sort(covp);
figure
subplot(1, 2, 1)
hold on
grid on
plot(nstimulus(infp:supp))
plot(nvalues(indp(1), infp:supp))
plot(nvalues(indp(end), infp:supp))
title('Stimulus, and neurons most and least correlated to positive stimulus', 'Interpreter', 'latex')
xlabel('Time index', 'Interpreter', 'latex')
ylabel('Normalized amplitude', 'Interpreter', 'latex')
legend('Stimulus', 'Neuron most correlated to negative stimulus', 'Neuron most correlated to positive stimulus')



%% Mean signal for two periods:

period = 200;
lper = length(nstimulus) / period;
stimean = mean(reshape(nstimulus, period, lper), 2);
signalp = mean(reshape(nvalues(indp(end), :), period, lper), 2);
signaln = mean(reshape(nvalues(indp(1), :), period, lper), 2);

% Plotting:
subplot(1, 2, 2)
hold on
grid on
plot(stimean)
plot(signaln)
plot(signalp)
title('Stimulus, and neurons most and least correlated to positive stimulus, AVERAGED', 'Interpreter', 'latex')
xlabel('Time index', 'Interpreter', 'latex')
ylabel('Mean normalized amplitude', 'Interpreter', 'latex')
legend('Stimulus', 'Neuron most correlated to negative stimulus', 'Neuron most correlated to positive stimulus')



%% Same analysis for differently correlated neurons:

lneu = length(covp);
interest = [0, 280, 805, 1920, 7200];
neuk = [1 + interest, ceil((lneu+1)/2), lneu - fliplr(interest)];

% Plotting:
figure
for i = 1:5
    % Whole signal:
    subplot(5, 2, 2*(i-1)+1)
    hold on
    grid on
    plot(nstimulus(infp:supp))
    plot(nvalues(indp(neuk(1+i)), infp:supp))
    if i ~= 5
        plot(nvalues(indp(neuk(end-i)), infp:supp))
    end
    title('Stimulus, and neurons most and least correlated to positive stimulus', 'Interpreter', 'latex')
    xlabel('Time index', 'Interpreter', 'latex')
    ylabel('Normalized amplitude', 'Interpreter', 'latex')
    % Averaged signal:
    signalp = mean(reshape(nvalues(neuk(1+i), :), period, lper), 2);
    signaln = mean(reshape(nvalues(neuk(end-i), :), period, lper), 2);
    subplot(5, 2, 2*i)
    hold on
    grid on
    plot(stimean)
    plot(signaln)
    plot(signalp)
    title('Stimulus, and neurons most and least correlated to positive stimulus, AVERAGED', 'Interpreter', 'latex')
    xlabel('Time index', 'Interpreter', 'latex')
    ylabel('Mean normalized amplitude', 'Interpreter', 'latex')
end



%% We are going to analyse the reaction of the neurons to the velocity & accelaration:

% Computing velocity and acceleration for stimulus:
nstimv = [0, nstimulus(2:end)-nstimulus(1:end-1)];
nstima = [0, nstimv(2:end)-nstimv(1:end-1)];
nvpos = [0, nstimpos(2:end)-nstimpos(1:end-1)];
napos = [0, nvpos(2:end)-nvpos(1:end-1)];
nvneg = [0, nstimneg(2:end)-nstimneg(1:end-1)];
naneg = [0, nvneg(2:end)-nvneg(1:end-1)];

% Computing correlations:
covposv = sum(nvalues.*nvpos, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nvpos.^2, 2));
covnegv = sum(nvalues.*nvneg, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(nvneg.^2, 2));
covpv = (covposv - mean(covposv)) ./ std(covposv);
covnv = (covnegv - mean(covnegv)) ./ std(covnegv);
covposa = sum(nvalues.*napos, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(napos.^2, 2));
covnega = sum(nvalues.*naneg, 2) ./ sqrt(sum(nvalues.^2, 2).*sum(naneg.^2, 2));
covpa = (covposa - mean(covposa)) ./ std(covposa);
covna = (covnega - mean(covnega)) ./ std(covnega);

% Plotting:
figure
subplot(1, 3, 1)
plot(covp, covn, '.', 'Color', [0, 0, 0.8])
grid on
axis equal
title('Correlations to positive \& negative position', 'Interpreter', 'latex')
xlabel('Correlation to positive position', 'Interpreter', 'latex')
ylabel('Correlation to negative position', 'Interpreter', 'latex')
subplot(1, 3, 2)
plot(covpv, covnv, '.', 'Color', [0.8, 0, 0])
grid on
axis equal
title('Correlations to positive \& negative velocity', 'Interpreter', 'latex')
xlabel('Correlation to positive velocity', 'Interpreter', 'latex')
ylabel('Correlation to negative velocity', 'Interpreter', 'latex')
subplot(1, 3, 3)
plot(covpa, covna, '.', 'Color', [0.8, 0.8, 0])
grid on
axis equal
title('Correlations to positive \& negative acceleration', 'Interpreter', 'latex')
xlabel('Correlation to positive acceleration', 'Interpreter', 'latex')
ylabel('Correlation to negative acceleration', 'Interpreter', 'latex')

% Whole plot:
totc = [covp, covn, covpv, covnv, covpa, covna];
colol = [[0, 0, 0.8]; [0, 0, 0.8]; [0.8, 0, 0]; [0.8, 0, 0]; [0.8, 0.8, 0]; [0.8, 0.8, 0]];
label = {'Positive position', 'Negative position', 'Positive velocity', 'Negative velocity', 'Positive acceleration', 'Negative acceleration'};

figure
for i = 1:6
    for j = i+1:6
        subplot(5, 5, 5*(i-1)+j-1)
        coltemp = (colol(i, :)+colol(j, :))/2;
        plot(totc(:, i), totc(:, j), '.', 'Color', coltemp)
        grid on
        axis equal
        if j == i+1
            ylabel(label{i}, 'Interpreter', 'latex')
            xlabel(label{i+1}, 'Interpreter', 'latex')
        end
    end
end



%% Trying to sum up a few neurons to obtain stability:

limitbl = 0.5;
biggest = (covp > limitbl);
lowest = (covn > limitbl);
figure
hold on
plot(nstimulus)
plot(nstimv)
toplot1 = mean(nvalues(biggest, :));
plot(toplot1+abs(min(toplot1)), 'Color', [0.5, 0.5, 0.5])
toplot2 = mean(nvalues(lowest, :));
plot(toplot2+abs(min(toplot2)), 'Color', [0, 0, 0])



%% Smoothing algorithm + peaks and comparison with speed:

% rank = 3;
% nsides = 12;
% addpath('~/Science/Hippolyte/rawDataAnalysis')
% toplot1 = toplot1 - min(toplot1);
% smooth1 = SGfilter(toplot1, rank, nsides);
% figure
% plot(toplot1, ':')
% hold on
% plot(smooth1)
% 
% peaks = peakDetection(smooth1, 'abs', 1);
% figure; 
% hold on
% plot(toplot1, ':')
% plot(smooth1)
% peaksf = zeros(3000, 1);
% peaksf(peaks) = max(smooth1);
% plot(peaksf)
% 
% vneur = [0, smooth1(2:end)-smooth1(1:end-1)];
% vneur = (vneur - mean(vneur)) / std(vneur);
% smoothv = SGfilter(vneur, rank, nsides);
% peaksv = peakDetection(smoothv, 'abs', 0.25);
% peaksfv = zeros(3000, 1);
% peaksfv(peaksv) = max(smoothv);
% figure
% hold on
% plot(vneur, ':')
% plot(smoothv)
% plot(peaksfv)
% plot(nstimulus)

nvelocity = [0, nstimulus(2:end)-nstimulus(1:end-1)];
nvelabs = abs(nvelocity);
nvelpos = nvelocity .* (nvelocity > 0);
nvelneg = -nvelocity .* (nvelocity < 0);
figure
subplot(5, 1, 1)
plot(nstimulus)
title('Normalized stimulus', 'Interpreter', 'latex')
grid on
subplot(5, 1, 2)
plot(nvelocity)
title('Velocity based on normalized stimulus', 'Interpreter', 'latex')
grid on
subplot(5, 1, 3)
plot(nvelabs)
title('Absolute value of velocity', 'Interpreter', 'latex')
grid on
subplot(5, 1, 4)
plot(nvelpos)
title('Positive velocity', 'Interpreter', 'latex')
grid on
subplot(5, 1, 5)
plot(nvelneg)
title('Negative velocity')
grid on

rank = 3;
nsides = 25;
addpath('~/Science/Hippolyte/rawDataAnalysis')
smoothvalues = SGfilter(nvalues, rank, nsides);
ndiffact = [zeros(size(smoothvalues, 1), 1), smoothvalues(:, 2:end)-smoothvalues(:, 1:end-1)];
ndiffact = (ndiffact - mean(ndiffact, 2)) ./ std(ndiffact, [], 2);
smoothpeed = SGfilter(ndiffact, rank, nsides);
figure
hold on
plot(nvalues(1, :), ':')
plot(smoothvalues(1, :))
plot(smoothpeed(1, :))


peakst = zeros(size(nvalues));
for i = 1:size(nvalues, 1)
    % First:
    peakstemp1 = peakDetection(smoothpeed(i, :), 'num', 30);
    peakstempf1 = zeros(1, 3000);
    peakstempf1(peakstemp1) = 1;
    % Second:
    peakstemp2 = peakDetection(smoothpeed(i, :), 'num', 60);
    peakstempf2 = zeros(1, 3000);
    peakstempf2(peakstemp2) = 1;
    % Third:
    peakstemp3 = peakDetection(smoothpeed(i, :), 'num', 120);
    peakstempf3 = zeros(1, 3000);
    peakstempf3(peakstemp3) = 0.5;
    peakst(i, :) = peakstempf1 + peakstempf2 + peakstempf3;
    if mod(i, 1000) == 0
        fprintf('Iteration %i out of %i \n', [i, size(nvalues, 1)])
    end
end
plot(peakst(1, :))
plot(nstimulus)
figure
image(peakst, 'CDataMapping', 'scaled')
colorbar

% Computing correlation:
cortot = zeros(size(nvalues, 1), 3);
cortot(:, 1) = sum(peakst.*nvelabs, 2) ./ sqrt(sum(peakst.^2, 2).*sum(nvelabs.^2, 2));
cortot(:, 2) = sum(peakst.*nvelpos, 2) ./ sqrt(sum(peakst.^2, 2).*sum(nvelpos.^2, 2));
cortot(:, 3) = sum(peakst.*nvelneg, 2) ./ sqrt(sum(peakst.^2, 2).*sum(nvelneg.^2, 2));
figure
image(cortot, 'CDataMapping', 'scaled')
colorbar
[~, indabs] = sort(cortot(:, 1));
figure
hold on
plot(peakst(indabs(end), :))
plot(nstimabs)




    









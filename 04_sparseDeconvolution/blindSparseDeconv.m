clear; close all; clc
% Matlab program for deconvolution paper.


%% Knowing parameters: 

% Constants:
a = 100;
b = 0.01;
dt = 0.01;
snr = 5;
sigma = a/snr;

% Spike train:
ltrain = 10*60;
fsig = 0.08;
fe = 2.5;
N = zeros(1, ltrain);
spikes = sin(fsig .* (dt:dt:ltrain*dt) .* 2 .* pi) + randn(1, ltrain) ;
% spikes = unique([1:40:ltrain, 3:41:ltrain, 2:42:ltrain, 4:25:ltrain]);
N(spikes > 1.5) = 1;

% Convolution kernel:
taud = 0.5;
taur = 0.05;
kern = @(t, taud, taur) (exp(-t./taud) - exp(-t./taur)) ./ max((exp(-(0:dt:3)./taud) - exp(-(0:dt:3)./taur)));

% Value of calcium signal:
F = zeros(size(N));
Fnonoise = zeros(size(N));
taunn = [];
for i = 1:length(N)
    Ntemp = N(1:i);
%     taunn = find(fliplr(Ntemp) == 1) - 1;
    if Ntemp(end) == 1
        taunn = [0, taunn+1];
    else
        taunn = taunn + 1;
    end
    intemp = sum(kern(taunn * dt, taud, taur)) * dt;
    F(i) = a * intemp + b + sigma * randn(1);
    Fnonoise(i) = a * intemp + b;
end
figure
hold on
t = dt:dt:ltrain*dt;
tc = dt/2:dt/2:ltrain*dt;
Nw = reshape([zeros(1, ltrain); N], 1, length(tc));
Nw = 0.1 * max(F) * Nw + b;
plot(tc, Nw, 'k')
plot(t, F, ':', 'Color', [1, 0.8, 0.2])
plot(dt:1/fe:ltrain*dt, F(1:1/(fe*dt):ltrain), 'o-', 'Color', [1, 0, 1])
grid on
legend('Spikes', 'Calcium signal', 'Acquired signal')
title('Spikes, calcium signal, and acquired signal', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

figure
hold on
plot(t, Fnonoise, 'Color', [0.8, 0.6, 0], 'LineWidth', 2)
plot(dt:1/fe:ltrain*dt, F(1:1/(fe*dt):ltrain), 'Color', [1, 0.7, 1])
Fsmooth = SGfilter(F(1:1/(fe*dt):ltrain), 5, 16);
plot(dt:1/fe:ltrain*dt, Fsmooth, 'k', 'LineWidth', 2)
title('Actual and acquired signal', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')


% %% Convolution matrix:
% sk = ltrain;
% K = zeros(sk);
% for i = 1:sk
%     for j = 1:sk
%         K(i, j) = dt * (i - j + 1);
%     end
% end
% K = kern(K);
% K = K .* (K > 0);
% Kinv = inv(K);
% figure
% subplot(1, 2, 1)
% image(K, 'CDataMapping', 'scaled')
% colorbar
% subplot(1, 2, 2)
% image(Kinv, 'CDataMapping', 'scaled')
% colorbar

%% Not knowing parameters:

% [dens, val] = hist(F);
% [~, madens] = max(dens);
% bg = val(madens);
% 
% Fneg = [F((F - bg) < 0), -F((F - bg) < 0)];
% sigmag = std(Fneg);



%% Analysing sinusoid:

% Comparison of std's:
cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-06-14/Run 03/Analysis/HDF5'
hippo = 'h5hippo.h5';
values = h5read(hippo, '/Data/Values');
stimulus = h5read(hippo, '/Data/Stimulus');

val1 = values(1, :);
med1 = median(val1); 
valmed1 = val1 - med1; 
valn1 = valmed1 .* (valmed1 < 0); 
sigma1 = std([valn1, -valn1]);

addpath('~/Science/Hippolyte/rawDataAnalysis')
out1 = SGfilter(valmed1, 3, 25);
diffv = valmed1 - out1;
sigmaout1 = std(diffv);
% 
% figure
% hold on
% plot(valmed1, ':')
% plot(out1)


% Sum of exponential kernels:
dt = 0.01;
t = 0:dt:5; taur = 0.1; taud = 0.5;
K = (exp(-t./taud) - exp(-t./taur)) ./ max(exp(-t./taud) - exp(-t./taur));
numb = 0:1000;
Ksum = K' .* numb;
figure
subplot(1, 2, 1)
image(Ksum, 'CDataMapping', 'scaled')
colorbar
subplot(1, 2, 2)
image(4 * (Ksum > 500) + 3 * (Ksum <= 500 & Ksum > 100) + 2 * (Ksum <= 100 & Ksum > 10) + (Ksum <= 10 & Ksum > 1), 'CDataMapping', 'scaled')
colorbar
% % Approximation:
% tsplit = ceil((taud*taur/(taud-taur)*log(taud/taur))/dt)*dt;
% Kapprox = zeros(1, length(t));
% for i = 1:length(t)
%     if t(i) < tsplit
%         Kapprox(i) = (1-exp(-t(i)/taur))/(1 - exp(-tsplit/taur));
%     else
%         Kapprox(i) = (1 - exp(-tsplit/taur) + exp(-t(i)/taud))/(1 - exp(-tsplit/taud) + exp(-tsplit/taur));
%     end
% end
% figure
% hold on
% plot(K, ':')
% plot(Kapprox)

% Getting exponentials:
val1 = SGfilter(F(1:1/(fe*dt):ltrain), 4, 22);
val1 = F(1:1/(fe*dt):ltrain);
taud = 0.5;
taur = 0.1;
texp = 0:0.4:max([5*taud, 0.4]);
kerndec = kern(texp, taud, taur);
spikes = zeros(1, length(val1));
for i = 1:length(val1)
    spikeskeep = spikes(max([1, i-length(kerndec)+2]):i);
    spikestemp = [zeros(1, length(kerndec)-length(spikeskeep) - 1), spikeskeep, 0];
    spikest = (val1(i) - sum(spikestemp .* fliplr(kerndec))) / kerndec(2);
    spikes(i) = spikest * (spikest > 0);
end
reconstitution = zeros(length(val1));
nkerndec = [kerndec(2:end), 0];
for i = 1:length(val1)
    spikestemp = [zeros(1, length(nkerndec)-i), spikes(max([1, i-length(nkerndec)+1]):i)];
    reconstitution(i) = sum(spikestemp .* fliplr(nkerndec));
end
figure
subplot(3, 1, 1:2)
hold on
plot(val1(1:500), ':')
plot(reconstitution(1:500))
subplot(3, 1, 3)
plot(spikes(1:500), 'k')






n = 0:45;
nx = 0:0.001:45;
lambda = [3, 9, 27]';
delta = 0.9;
espe = lambda .* delta;
expo = exp(-nx .* 1 ./ espe) ./ espe;
poisson = exp(-espe) .* (espe .^ n) ./ factorial(n);
figure
subplot(2, 1, 1)
hold on
for i = 1:length(lambda)
    plot(n, poisson(i, :), 'o:')
end
subplot(2, 1, 2)
hold on
for i = 1:length(lambda)
    plot(nx, expo(i, :))
end


%% Presentation: 

% Constants:
a = 100;
b = 0.01;
dt = 0.01;
snr = 50;
sigma = a/snr;

% Spike train:
ltrain = 10*60;
fsig = 0.2;
fe = 2.5;
N = zeros(1, ltrain);
spikes = sin(fsig .* (dt:dt:ltrain*dt) .* 2 .* pi) + randn(1, ltrain) ;
% spikes = unique([1:40:ltrain, 3:41:ltrain, 2:42:ltrain, 4:25:ltrain]);
N(spikes > 1.5) = 1;

% Convolution kernel:
taud = 0.5;
taur = 0.05;
kern = @(t, taud, taur) (exp(-t./taud) - exp(-t./taur)) ./ max((exp(-(0:dt:3)./taud) - exp(-(0:dt:3)./taur)));

% Value of calcium signal:
F = zeros(size(N));
Fnonoise = zeros(size(N));
taunn = [];
for i = 1:length(N)
    Ntemp = N(1:i);
%     taunn = find(fliplr(Ntemp) == 1) - 1;
    if Ntemp(end) == 1
        taunn = [0, taunn+1];
    else
        taunn = taunn + 1;
    end
    intemp = sum(kern(taunn * dt, taud, taur)) * dt;
    F(i) = a * intemp + b + sigma * randn(1);
    Fnonoise(i) = a * intemp + b;
end
figure
t = dt:dt:ltrain*dt;
tc = dt/2:dt/2:ltrain*dt;
Nw = reshape([zeros(1, ltrain); N], 1, length(tc));
Nw = 0.1 * max(F) * Nw + b;
subplot(4, 1, 1)
plot(tc, Nw, 'k')
grid on
title('Spike train', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')
subplot(4, 1, 2:4)
plot(t, F, 'Color', [1, 0.8, 0.2])
grid on
title('Calcium signal', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')

figure
subplot(4, 1, 1:3)
plot(t, F, 'Color', [1, 0.8, 0.2])
grid on
title('Calcium signal', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')
subplot(4, 1, 4)
plot(tc, Nw, 'k')
grid on
title('Spike train', 'Interpreter', 'latex')
xlabel('Time [s]', 'Interpreter', 'latex')


t = 0:0.4:2.8;
gg = kern(t, taud, taur);
for i = 1:30; K(i:i+7, i) = gg'; end
K = K(1:30, :);
N = rand(1, 30);
N(N > 0.7) = 1;
N(N <= 0.7) = 0;
t = 0:0.4:29*0.4;
figure; 
subplot(4, 4, 1:3); Np = reshape([N; zeros(size(N))], 1, 60); tp = 0:0.2:59*0.2; plot(tp, Np, 'k'); 
axis([0, 59*0.2+0.8, 0, 1.5])
subplot(4, 4, [5, 6, 7, 9, 10, 11, 13, 14, 15]); image(K, 'CDataMapping', 'scaled'); colorbar; 
subplot(4, 4, [8, 12, 16]); plot(K*N', t, 'r'); set(gca,'ydir','reverse'); grid on
axis([0, max(K*N')+0.1, 0, max(t)])



%% Trying genSetConv:

params.nex = 10;
params.pts = 500;
params.fs = 2.5;
params.fr = 1000;
params.taur = [0.01, 0.1];
params.taud = [0.1, 1];
params.a = 100;
params.b = -10;
params.noise = 10;
params.ssec = 100;
params.per = 30;
params.perinf = 0.01;
params.indic = 10;
[data, info, N] = genSetConv(params);
figure
hold on
for i = 1:size(data, 2)
    plot(0:params.pts-1, data(:, i))
end



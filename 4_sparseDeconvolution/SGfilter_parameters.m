clear; close all; clc
warning('off')
tic


%% Knowing parameters: 

% Constants:
a = 100;
b = 0.01;
dt = 0.01;
minute = 1;
ntest = 40;

% Parameters:
snr = [1, 2.5, 5, 10, 25, 50, 100];
fsig = (0.005:0.005:0.1);
fe = [0.25, 0.5, 1, 2.5, 5];
param1 = 7;
param2 = 30;
% snr = [1, 100];
% fsig = [0.05, 1];
% fe = [0.1, 10];

% Final matrices:
lsnr = length(snr);
lfsig = length(fsig);
lfe = length(fe);
indpol = zeros(lsnr, lfsig, lfe, ntest);
numpt = indpol;
minstep = indpol;
minsmooth = zeros(param1, param2);

% Spike train:
ltrain = 60*minute/dt;
N = zeros(1, ltrain);
t = dt:dt:ltrain*dt;

% Convolution kernel:
taud = 0.5;
taur = 0.1;
kern = @(t, taud, taur) (exp(-t./taud) - exp(-t./taur)) ./ max((exp(-(0:dt:3)./taud) - exp(-(0:dt:3)./taur)));


% Main algorithm:
for isnr = 1:lsnr
    % Noise:
    sigma = a ./ snr(isnr);
    for ifsig = 1:lfsig
        for ife = 1:lfe
            te = 1:1/(fe(ife)*dt):ltrain;
            for itest = 1:ntest
                % Spike train:
                spikes = sin(fsig(ifsig) .* t .* 2 .* pi) + randn(1, ltrain) ;
                N(spikes > 1.5) = 1;
                % Value of calcium signal:
                Fnonoise = zeros(size(N));
                F = Fnonoise;
                taunn = [];
                for in = 1:length(N)
                    Ntemp = N(1:in);
                    if Ntemp(end) == 1
                        taunn = [0, taunn+1];
                    else
                        taunn = taunn + 1;
                    end
                    intemp = sum(kern(taunn * dt, taud, taur)) * dt;
                    Fnonoise(in) = a * intemp + b;
                    F(in) = a * intemp + b + sigma * randn(1);
                end
                % Optimization:
                for k1 = 1:param1
                    param2bis = floor((length(te)-1)/2);
                    for k2 = 1:min([param2, param2bis])
                        Fsmooth = SGfilter(F(te), k1-1, k2);
                        minsmooth(k1, k2) = sum((Fnonoise(te)-Fsmooth).^2);
                    end
                end
                minsmooth(minsmooth == 0) = max(minsmooth(:));
                [minstep(isnr, ifsig, ife, itest), indmin] = min(minsmooth(:));
                [indpol(isnr, ifsig, ife, itest), numpt(isnr, ifsig, ife, itest)] = ind2sub(size(minsmooth), indmin);
            end
            fprintf('%d %% done in %d seconds \n', [100*(ntest*lfe*lfsig*(isnr-1) + ntest*lfe*(ifsig-1) + ntest*ife)/(lsnr*lfsig*lfe*ntest), toc]);
        end
    end
end




nummean = mean(numpt, 4);
numstd = std(numpt, [], 4);
indmean = mean(indpol, 4);
indstd = std(indpol, [], 4);

nummean2d = zeros(lsnr, lfsig*lfe);
indmean2d = zeros(lsnr, lfsig*lfe);
for i = 1:lfe
    nummean2d(:, (lfsig*(i-1)+1):(lfsig*i)) = nummean(:, :, i);
    indmean2d(:, (lfsig*(i-1)+1):(lfsig*i)) = indmean(:, :, i);
end

figure
subplot(2, 1, 1)
hold on
image(indmean2d, 'CDataMapping', 'scaled')
for i = 1:lfe-1
    line([i*lfsig+0.5, i*lfsig+0.5], [0.5, lsnr+0.5], 'Color', [0, 0, 0])
end
colorbar
subplot(2, 1, 2)
image(nummean2d, 'CDataMapping', 'scaled')
for i = 1:lfe-1
    line([i*lfsig+0.5, i*lfsig+0.5], [0.5, lsnr+0.5], 'Color', [0, 0, 0])
end
colorbar

load('indpol.mat');
load('numpt.mat');
SGparams.snr = snr;
SGparams.fsig = fsig;
SGparams.fe = fe;
SGparams.ntest = ntest;
SGparams.indpol = indpol;
SGparams.numpt = numpt;







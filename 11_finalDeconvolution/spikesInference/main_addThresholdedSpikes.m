clear; close all; clc
%% This code generates spikes and thresholded spikes for an HDF5 file. %%

% -ENTER IN THIS SECTION path to BSD folder.
%
% -ENTER PATHS to HDF5 in the next section, as a cell (ptemp).
% -ENTER PATH to taur_taud.txt file in the next session, as a char (psave)
%    -If file does not exist, enter path to the folder that will get file.
%    -If you created this file outside of this program, use the following
%     convention: 'tau_rise {{tau_rise}} tau_decay {{tau_decay}}'.
% -IF NEEDED to compute the DFFaligned, change value two sections below 
%  (do_aligned) to true, if needed to add it to the HDF5, change value two 
%  sections below (add_aligned) to true.
%
% -YOU CAN CHANGE parameters in estimateTimeConstantsHDF5 three sections
%  below. Line 18 of the function show possible inputs. You can also change
%  these parameters four sections below.
%
% -IF LABELS ARE MISSING from HDF5 file, and you want an estimation based
%  on particular labels, please enter path to MaskDatabase.mat on line 7 
%  of addLabels.m. MaskDatabase.mat can be downloaded at this address:
%  https://engertlab.fas.harvard.edu/Z-Brain/download
%
% -SIDE NOTE: This code uses the HDF5 architecture defined by laboratoire 
%  Jean Perrin.

% Path to BSD algorithm:
addpath(genpath('/home/ljp/Programs/BSD'))



%% Provide HDF5 path:

% Provide as a cell all the paths to HDF5 necessary to compute tau rise and 
% tau decay. Inferred spikes will be computed for all these HDF5:
ptemp = {'/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2019-03-26/Run 10/Analysis/HDF5/2019-03-26(Run10).h5';
          '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2019-03-26/Run 09/Analysis/HDF5/2019-03-26(Run09).h5'};
% Provide path to which taur_taud.txt file will be saved:
psave = '/home/ljp/Science/Hippolyte/Small-programs/11_finalDeconvolution/spikesInference';
      


%% Interpolate to find DF_aligned:

% Compute DFFaligned? If true proceding:
do_aligned = false; % CHANGE IF NEEDED
if do_aligned
    for p = 1:length(ptemp)
        times = h5read(ptemp{p}, '/Data/Brain/Times');
        delays = h5read(ptemp{p}, '/Data/Brain/TimeDelays');
        dff = h5read(ptemp{p}, '/Data/Brain/Analysis/DFF');
        translaTime = times + delays;
        DFFaligned = zeros(size(dff));
        for i = 1:size(dff, 1)
            DFFaligned(i, :) = interp1(translaTime(i, :), dff(i, :), times);
            showProgress(i, size(dff, 1)); % DEACTIVATE IF NEEDED
        end
        fprintf('\n');
        % Add DFFaligned to HDF5? If true proceding:
        add_aligned = false; % CHANGE IF NEEDED
        if add_aligned
            h5create(ptemp{p}, '/Data/Brain/Analysis/DFFaligned', size(DFFaligned), 'Datatype', 'single');
            h5write(ptemp{p}, '/Data/Brain/Analysis/DFFaligned', cast(DFFaligned, 'single'))
        end
    end
end



%% Program defines constants based on tectum region:

% Check if constants already exist:
taufile = fullfile(psave, 'taur_taud.txt');
if ~exist(taufile, 'file')
    % Define constants:
    fprintf('Launching algorithm. \n');
    tau = cell(length(ptemp), 2);
    for p = 1:length(ptemp)
        [taur, taud] = estimateTimeConstantsHDF5(ptemp{p}, 'labels', 106, 'parfor', true, 'aligned', true); % CHANGE PARAMETERS IF NEEDED
        tau{p, 1} = taur;
        tau{p, 2} = taud;
    end
    % Compute final constants:
    taur = 0;
    taud = 0;
    for p = 1:length(ptemp)
        taur = taur + (tau{p, 1}/length(ptemp));
        taud = taud + (tau{p, 2}/length(ptemp));
    end
    % Save as text file:
    fid = fopen(taufile, 'wt');
    fprintf(fid, 'tau_rise %f tau_decay %f', [taur, taud]);
else
    % Else recover values:
    fid = fopen(taufile, 'r');
    tau_temp = fscanf(fid, '%c');
    tau_temp = split(string(tau_temp));
    taur = double(tau_temp(2));
    taud = double(tau_temp(4));
end



%% Computing inferred and thresholded spikes:

for p = 1:length(ptemp)
    % BSD parameters in random order (refer to estimateTimeConstantsHDF5.m):
    Palg = struct;
    Palg.tauRise = taur;
    Palg.tauDecay = taud;
    Oalg = struct;
    Oalg.adaptive = 1; 
    Oalg.iterations = 5; 
    dff = h5read(ptemp{p}, '/Data/Brain/Analysis/DFFaligned'); % can also put '/Data/Brain/Analysis/DFF'  
    dff = dff';  
    dff = cast(dff, 'double'); 
    time = h5read(ptemp{p}, '/Data/Brain/Times'); 
    time = cast(time, 'double'); 
    dff(abs(dff) > 10*std(dff(:))) = 1e-3; 
    dff(isnan(dff)) = 1e-3; 
    dff(dff == 0) = 1e-3; 
    Oalg.Time = size(dff, 1);
    Oalg.dt = mean(gradient(time));
    Oalg.nNeurons = size(dff, 2);
    % Go through BSD:
    fprintf('Launching BSD algorithm... \n');
    [N, ~, ~, Pphys] = pBSD(dff, Oalg, Palg); % can also use simple BSD (comment next line then)
    delete(gcp('nocreate'));
    % Threshold:
    N = N';
    thresh = Pphys.threshold';
    Nthresh = (N > thresh) * 1;
    % Add to HDF5:
    % Inferred spikes:
    h5create(ptemp{p}, '/Data/Brain/Analysis/InferredSpikes', size(N), 'Datatype', 'single');
    h5write(ptemp{p}, '/Data/Brain/Analysis/InferredSpikes', cast(N, 'single'));
    h5writeatt(ptemp{p}, '/Data/Brain/Analysis/InferredSpikes', 'tauRise', taur);
    h5writeatt(ptemp{p}, '/Data/Brain/Analysis/InferredSpikes', 'tauDecay', taud);
    % Thresholded spikes:
    h5create(ptemp{p}, '/Data/Brain/Analysis/ThresholdedSpikes', size(Nthresh), 'Datatype', 'single');
    h5write(ptemp{p}, '/Data/Brain/Analysis/ThresholdedSpikes', cast(Nthresh, 'single'));
    h5writeatt(ptemp{p}, '/Data/Brain/Analysis/ThresholdedSpikes', 'tauDecay', taud);
    h5writeatt(ptemp{p}, '/Data/Brain/Analysis/ThresholdedSpikes', 'tauRise', taur);
end



clear; close all; clc
% Path to BSD algorithm:
addpath(genpath('/home/ljp/Programs'))
% Path to function:
addpath(genpath('/home/ljp/Science/Hippolyte/multiSensorint'))
% Path to showProgress:
addpath(genpath('/home/ljp/Science/Hippolyte/Small-programs'))
% This code generates spikes and thresholded spikes for an HDF5 file.



%% Provide HDF5 path:

ptemp = '/home/ljp/Science/GeoffreysComputer/Projects/RLS/Data/2019-03-26/Run 09/Analysis/HDF5/2019-03-26(Run09).h5';
ptemp = '/home/ljp/Science/Hippolyte/ALL_DATASETS/2018-05-24Run07.h5';



%% Interpolate to find DF_aligned:

% Compute DFFaligned? If true proceding:
do_aligned = true;
if do_aligned
    times = h5read(ptemp, '/Data/Brain/Times');
    delays = h5read(ptemp, '/Data/Brain/TimeDelays');
    dff = h5read(ptemp, '/Data/Brain/Analysis/DFF');
    translaTime = times + delays;
    DFFaligned = zeros(size(dff));
    for i = 1:size(dff, 1)
        DFFaligned(i, :) = interp1(translaTime(i, :), dff(i, :), times);
%         showProgress(i, size(dff, 1));
        showProgress_sign(i, size(dff, 1), 50, 50, '>>>', ' ');
%         showProgress_signCheck(i, size(dff, 1), 'progress_sign', '>>>', 'progress_bar', ' ');
    end
    fprintf('\n');
    % Add DFFaligned to HDF5? If true proceding:
    add_aligned = false;
    if add_aligned
        h5create(ptemp, '/Data/Brain/Analysis/DFFaligned', size(DFFaligned), 'Datatype', 'single');
        h5write(ptemp, '/Data/Brain/Analysis/DFFaligned', cast(DFFaligned, 'single'))
    end
end



%% Program defines constants based on tectum region:

fprintf('Launching algorithm. \n');
[taur, taud] = estimate_time_constants_HDF5(ptemp, 'labels', 106, 'parfor', false, "aligned", true);



%% Taking 

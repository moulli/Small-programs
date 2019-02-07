%% TEMPLATE : do not change this file !!!

% if you want to change this file
% create instead a copy in other file (ex TEMPLATE_Hugo.m)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
%                                                                       %
%                         DO NOT EDIT ANYMORE !                         %
%                                                                       %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear; clc

%% add programs

addPrograms('/home/ljp/');

%% sample focus

% root = '/home/ljp/Science/Hippolyte/dataHugo/';
% study = '';
% date = '2018-11-21';
% run = 'Run 03';
root = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/';
study = '';
date = '2018-05-24';
run = 'Run 07';


F = NT.Focus(root, study, date, run);

% sample parameters

clear Analysis
F
Analysis.Layers = 2:10;             % Layers to analyse
Analysis.RefIndex = 1;              % index of the reference frame for drift correction
Analysis.RefStack = '';             % external reference stack if exists
Analysis.BaselineWindow = 50;       % time in seconds of the baseline window
Analysis.BaselinePercentile = 10;   % percentile for baseline computation
Analysis.Lineage = 'Nuclear';       % possible values : 'Nuclear', 'Cytoplasmic'
Analysis.StimulusFrequency = 0.2;   % frequency of stimulus (Hz) for phasemap computation
Analysis.Stimulus = 'sinus';        % type of stimulus (step/sinus)
Analysis.Overwrite = false;         % true if program is allowed to overwrite
Analysis.RefBrain = 'zBrain_Elavl3-H2BRFP_198layers.nhdr'; % choose refbrain to map onto

% loads the parameters in the current focus
F.Analysis = Analysis;

%% sample viewer (collection of all viewer functions)

Focused.stackViewer(F, 'source');                   % view raw stack
Focused.stackViewer(F, 'ROImask');                  % view the user defined ROI
seeDriftCorrection(F, 5);                           % plays a movie of corrected images for layer 5                                  
Focused.stackViewer(F, 'corrected', [400 1200]);    % plots the corrected stack with values between 400 and 1200
Focused.stackViewer(F, 'graystack');                % view the graystack
stackViewer2D(F, 'BaselineNeuron');                 % view the baseline computed per neuron
stackViewer2D(F, 'DFFNeuron');                      % view the dff computed per neuron
stackViewer2D(F, 'BaselinePixel');                  % view the baseline computed per pixel
stackViewer2D(F, 'DFFPixel');                       % view the dff computed per pixel
Focused.phaseMapViewer(F, 'signal pixel',1000);     % view the phasemap computed on signal per pixel
Focused.phaseMapViewer(F, 'dff pixel',10);          % view the phasemap computed on dff per pixel
Focused.phaseMapViewer(F, 'dff neuron');            % view the phasemap computed on dff per neuron

%% sample functions (to run the analysis function by function)

semiAutoROI(F);                     % lets the user select the ROI
Focused.driftCompute(F, 'both')     % computes the fast drift then the slow drift
driftApply(F);                      % applies computed drift
driftComputeAndApply(F, 'on');      % calculates the drift for every layer independently

computeBackground(F);               % computes the background value for each layer
createGrayStack(F);                 % creates a gray stack by averaging on image each 77 
segmentBrain(F, 'graystack','RC');  % segment brain using the Raphaël Candelier function

% --- per neuron
computeBaseline(F, 'neuron');       % computes baseline per neuron
computeDFF(F, 'neuron');            % computes dff per neuron
computePhaseMap(F, 'neuron');       % computes phasemap per neuron

% --- per pixel
computeBaseline(F, 'pixel');        % computes baseline per pixel
computeDFF(F, 'pixel');             % computes dff per pixel
phaseMapPixel(F);                   % computes phasemap per pixel on dff
computePhaseMap(F, 'pixel', 'signal');  % computes baseline per pixel on signal
computePhaseMap(F, 'pixel', 'dff');     % computes baseline per pixel on dff (Geoffrey's function)

% --- export
mapToRefBrain(F, 'affine', '', 'graystack');                % maps graystack to reference brain (affine)
mapToRefBrain(F, 'convertcoord', 'affine', 'graystack');    % converts coordinates of segmented neurons
exportToHDF5(F);               % export data ub hdf5 file









%% Analysis Hippolyte:

% Creating mmap for corrected drift:
m = Focused.Mmap(F, 'corrected');

% Computing std:
tstep = 1;
mxinf = 14;
mxlen = m.x-mxinf+1;
myinf = 6;
mysup = m.y-10;
mylen = mysup-myinf+1;
mcor = zeros(mxlen, mylen, m.z);
        
stot = mylen * m.z;    
for j = myinf:mysup
    for k = 1:m.z
        timetemp = permute(m(mxinf:m.x, j, m.Z(k), tstep:tstep:m.t), [4, 1, 2, 3]);
        timestd = std(double(timetemp));
        mcor(:, j-myinf+1, k) = permute(timestd, [2, 1, 3]);
        count = (j-1)*m.z+k;
        if mod(count, ceil(stot/10)) == 0
            fprintf('%2.1f %% done \n', count/ceil(stot/100));
        end
    end
end

couche = 5;
figure
image(mcor(:, :, couche), 'CDataMapping', 'scaled')
title('STD matrix from raw data', 'Interpreter', 'latex')
axis equal
colorbar

mstdf = zeros(size(mcor));
division = [10, 15, 20, 50];
for i = 1:length(division)
    mstdf(mcor > division(i)) = i;
end
figure
image(mstdf(:, :, couche), 'CDataMapping', 'scaled')
title('Steps of STD matrix from raw data', 'Interpreter', 'latex')
axis equal
colorbar

% Computing 3d std of mstd to highligh zones of high change:
[msi, msj, msk] = size(mcor);
mstdg3 = zeros(size(mcor));
for i = 1:msi
    for j = 1:msj
        for k = 1:msk
            neigh = [i + [-1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1], ...
                     j + [-1; -1; -1; 0; 0; 0; 1; 1; 1; -1; -1; -1; 0; 0; 0; 1; 1; 1; -1; -1; -1; 0; 0; 0; 1; 1; 1], ...
                     k + [-1; -1; -1; -1; -1; -1; -1; -1; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1]];
            neigh = neigh(all(neigh, 2) & neigh(:, 1) <= msi & neigh(:, 2) <= msj & neigh(:, 3) <= msk, :);
            neighind = sub2ind(size(mcor), neigh(:, 1), neigh(:, 2), neigh(:, 3));
            mstdg3(i, j, k) = std(mcor(neighind));
            count = (i-1)*mj*mk+(j-1)*mk+k;
            if mod(count, ceil(mi*mj*mk/100)) == 0
                fprintf('%2.1f %% done \n', count/ceil(mi*mj*mk/100));
            end
        end
    end
end
figure
image(mstdg3(:, :, couche), 'CDataMapping', 'scaled')
title('STD of difference from STD matrix for raw data (3D)', 'Interpreter', 'latex')
axis equal
colorbar

% mstdN = (mcor - mean(mcor(:))) ./ std(mcor(:));
% mstdg3N = (mstdg3 - mean(mstdg3(:))) ./ std(mstdg3(:));
% sqdmstd = (mstdN - mstdg3N) .^ 2;
% figure
% image(sqdmstd(300:500, 200:400, couche), 'CDataMapping', 'scaled')
% title('Difference between STD matrix and STD difference from STD matrix (3D)', 'Interpreter', 'latex')
% axis equal
% colorbar

% Computing 2d std of mstd to highligh zones of high change:
[msi, msj, msk] = size(mcor);
mstdg2 = zeros(size(mcor));
for i = 1:msi
    for j = 1:msj
        for k = 1:msk
            neigh = [i + [-1; 0; 1; -1; 0; 1; -1; 0; 1], ...
                     j + [-1; -1; -1; 0; 0; 0; 1; 1; 1], ...
                     k + [0; 0; 0; 0; 0; 0; 0; 0; 0]];
            neigh = neigh(all(neigh, 2) & neigh(:, 1) <= msi & neigh(:, 2) <= msj & neigh(:, 3) <= msk, :);
            neighind = sub2ind(size(mcor), neigh(:, 1), neigh(:, 2), neigh(:, 3));
            mstdg2(i, j, k) = std(mcor(neighind));
            count = (i-1)*msj*msk+(j-1)*msk+k;
            if mod(count, ceil(msi*msj*msk/100)) == 0
                fprintf('%2.1f %% done \n', count/ceil(msi*msj*msk/100));
            end
        end
    end
end
figure
image(mstdg2(:, :, couche), 'CDataMapping', 'scaled')
title('STD of difference from STD matrix for raw data (2D)', 'Interpreter', 'latex')
axis equal
colorbar


% Using mstdg2 as a mask (bof):
mstdg2N = (mstdg2 - mean(mstdg2(:))) ./ std(mstdg2(:));
mstdg2m = mstdg2N;
mstdg2m(mstdg2m <= 0) = 1000 * max(mstdg2m(:));
mstdg2m = mstdg2m + 1;
mstdMasked = mstdN ./ mstdg2m;
figure
image(mstdMasked(:, :, couche), 'CDataMapping', 'scaled')
title('Mask for neurons based on std', 'Interpreter', 'latex')
axis equal
colorbar

% Getting rid of background on original data:
mstdg2m2 = mstdg2;
mstdg2m2(mstdg2m2 <= 0.5) = 0;
mstdMasked2 = mstd .* (mstdg2m2 > 0.5);
figure
image(mstdMasked2(:, :, couche), 'CDataMapping', 'scaled')
title('Original data with mask getting rid of background', 'Interpreter', 'latex')
axis equal
colorbar


% % Correlation from m:
% tstep = 1000;
% mxinf = 14;
% mxsup = m.x;
% mxlen = mxsup-mxinf+1;
% myinf = 6;
% mysup = m.y-10;
% mylen = mysup-myinf+1;
% mcor = zeros(mxlen, mylen, m.z);
% 
% stot = mxlen * mylen * m.z;   
% tic
% for j = myinf:mysup
%     for k = 1:m.z
%         kz = m.Z(k);
%         yzkeep = [j + [-1; 0; 1; -1; 1; -1; 0; 1];
%                   kz + [-1; -1; -1; 0; 0; 0; 1; 1; 1]];
%         neigh = neigh(all(neigh, 2) & neigh(:, 1) <= mysup & neigh(:, 2) <= max(m.Z), :);
%         mtemp = double(m(mxinf:mxsup, min(neigh(:, 1)):max(neigh(:, 1)), min(neigh(:, 2)):max(neigh(:, 2))));
%         for i = mxinf:mxsup
%             neigh2 = [i + [-1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1], ...
%                       j + [-1; -1; -1; 0; 0; 0; 1; 1; 1; -1; -1; -1; 0; 0; 1; 1; 1; -1; -1; -1; 0; 0; 0; 1; 1; 1], ...
%                       kz + [-1; -1; -1; -1; -1; -1; -1; -1; -1; 0; 0; 0; 0; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1]];
%             neigh2 = neigh2(all(neigh2, 2) & neigh2(:, 1) <= mxsup & neigh2(:, 2) <= mysup & neigh2(:, 3) <= max(m.Z), :);
%             neigh2(:, 3) = neigh2(:, 3) - min(neigh2(:, 3)) + min(m.Z);
%             mneigh = zeros(size(neigh2, 1), size(m(i, j, kz, :), 4));
%             for g = 1:size(neigh2)
%                 mneigh(i, :) = permute(double(m(neigh2(g, 1), neigh2(g, 2), neigh2(g, 3), :)), [1, 4, 2, 3]);
%             end
%             ownm = permute(double(m(i, j, kz, :)), [1, 4, 2, 3]);
%             mcor(i-mxinf+1, j-myinf+1, k) = mean(sum(ownm.*mneigh, 2) ./ sqrt(sum(ownm.^2, 2).*sum(mneigh.^2, 2)));
%             count = (j-myinf)*mxlen*m.z+(k-1)*mxlen+i-mxinf+1;
%             if mod(count, ceil(stot/10000)) == 0
%                 fprintf('%2.1f %% done, in %.0f seconds \n', [count/ceil(stot/10000), toc]);
%             end
%         end
%     end
% end










%%






% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
%                                                                       %
%                         DO NOT EDIT ANYMORE !                         %
%                                                                       %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %








%% adding programs

function addPrograms(root)
%addPrograms adds the matlab programs to the path and loads the caTools library
% root is the root of the programs (ex /home/ljp/programs)

    % adds matlab programs path
    addpath(genpath(fullfile(root,'Programs', 'easyRLS','Matlab')));
    addpath(genpath(fullfile(root,'Programs', 'NeuroTools','Matlab')));

    dir = NT.Focus.architecture(root, 'none');

    if ismac
        warning('test if caTools is ok for mac');
    elseif isunix
        cd(dir('caTools'))
        [~,~] = loadlibrary('caTools.so',...
                            'caTools.h');
    elseif ispc
        cd(dir('caTools'))
        [~,~] = loadlibrary('caTools.dll',...
                            'caTools.h');
    else
        disp('Platform not supported')
    end
    
    cd(root);
    disp('done');
end
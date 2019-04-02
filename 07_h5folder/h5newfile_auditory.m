clear; close all; clc
addpath(genpath('/home/ljp/Programs/'));
addpath(genpath('/home/ljp/Science/Hippolyte'));



%% Defining path:

Fpath = '/home/ljp/Science/Hippolyte/auditory_h5/auditory_20150303Run07.h5';
times = h5read(Fpath, '/Data/Times');
stimulus = h5read(Fpath, '/Data/stimulus');
figure; plot(times(1:100), stimulus(1:100))



%% Trying function h5switch:

nfolder = struct;
nfolder.path = '/home/ljp/Science/Hippolyte/auditory_h5';
nfolder.name = 'newHDF5_auditory';
nfolder.filename = '2015-03-03Run07_a.h5';

information = struct;
% information.description = 'test to see.';
information.autodescription = true;
information.date = '2015-03-03';
information.run = 7;
% information.rate = 3;
% information.layers = 17;
% information.increment = 4; %in micrometers
% information.line = 'poisson random';
% information.age = 7;
% information.autoID = true;
stim.name = {'auditory1'};
stim.sensorytype = {'auditory'};
stim.stimtype = {'pulses'};
stim.frequency = {0.1};
information.stimulus = stim;
behav.name = {};
behav.type = {};
behav.frequency = {};
information.behaviour = behav;

overwrite = struct;
overwrite.folder = false;
overwrite.file = true;

h5switch_auditory(Fpath, nfolder, information, overwrite);



%% Opening new HDF5 file:

% hdf5test = fullfile(nfolder.path, nfolder.name, nfolder.filename);
% h5disp(hdf5test)
% 
% dff = h5read(hdf5test, '/Data/Brain/Analysis/DFF');
% stim = h5read(hdf5test, '/Data/Stimulus/vestibular1/motorAngle');
% dffanim = 0;
% if dffanim
%     figure('units','normalized','outerposition',[0 0 1 1])
% %     [~, indff] = sort(var(dff, [], 2), 'descend');
% %     dffsort = dff(indff, :);
%     mdff = mean(dff);
%     for i = 1:3000
%         pause(0.000001)
%         subplot(5, 1, 1:3)
%         plot(dff(:, i), '.')
%         title(['Iteration ', num2str(i)])
%         axis([1, size(dff, 1), -0.2, 1])
%         subplot(5, 1, 4)
%         hold off
%         plot(mdff)
%         hold on
%         plot(i, mdff(i), 'o', 'MarkerFaceColor', 'r')
%         subplot(5, 1, 5)
%         hold off
%         plot(stim, 'y')
%         hold on
%         plot(i, stim(i), 'or', 'MarkerFaceColor', 'r')
%     end
% end



%% Guillaume's data:

% addpath('/home/ljp/Science/Guillaume/Thermotaxis/Datasets')
% h5disp('20171115_Run06_rp_Tset=27.h5')
% 
% 
% 
% %% Working on fillMetaHDF5:
% 
% % filename = fullfile(F.dir('HDF5'), [erase(F.name, ' ') '.h5']);
% % filename = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run07';
% addpath(genpath('/home/ljp/Programs'))
% information = fillMetaHDF5(F);
% % 
% % attstruct = h5info(filename, '/Metadata')
% % attstruct.Attributes.Name









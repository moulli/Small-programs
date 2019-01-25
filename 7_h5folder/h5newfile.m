clear; close all; clc



%% Defining focus:

addpath(genpath('/home/ljp/Programs/'));

root = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/';
study = '';
date = '2018-05-24';
run = 'Run 08';

F = NT.Focus(root, study, date, run);



%% Trying function h5switch:

nfolder = struct;
nfolder.path = '/home/ljp/Science/Hippolyte';
nfolder.name = 'newHDF5';
nfolder.filename = 'testfile1.h5';

information = struct;
% information.description = 'test to see.';
information.autodescription = true;
information.date = '2018-05-24';
information.run = 6;
information.rate = 3;
information.layers = 17;
information.increment = 4; %in micrometers
information.line = 'poisson random';
information.age = 7;
information.autoID = true;
stim.name = {'vestibular1', 'thermotaxis1'};
stim.sensorytype = {'vestibular', 'thermotaxis'};
stim.stimtype = {'sinus', 'step'};
stim.frequency = {20, 0.4};
information.stimulus = stim;
behav.name = {'eyesmovement1', 'tailmovement1'};
behav.type = {'eye tracking', 'tail tracking'};
behav.frequency = {100, 80};
information.behaviour = behav;

overwrite = struct;
overwrite.folder = false;
overwrite.file = false;

h5switch(F, nfolder, information, overwrite);



%% Opening new HDF5 file:

hdf5test = fullfile(nfolder.path, nfolder.name, nfolder.filename);
h5disp(hdf5test)

dff = h5read(hdf5test, '/Data/Brain/Analysis/DFF');
stim = h5read(hdf5test, '/Data/Stimulus/vestibular1/motorAngle');
dffanim = 0;
if dffanim
    figure('units','normalized','outerposition',[0 0 1 1])
%     [~, indff] = sort(var(dff, [], 2), 'descend');
%     dffsort = dff(indff, :);
    mdff = mean(dff);
    for i = 1:3000
        pause(0.000001)
        subplot(5, 1, 1:3)
        plot(dff(:, i), '.')
        title(['Iteration ', num2str(i)])
        axis([1, size(dff, 1), -0.2, 1])
        subplot(5, 1, 4)
        hold off
        plot(mdff)
        hold on
        plot(i, mdff(i), 'o', 'MarkerFaceColor', 'r')
        subplot(5, 1, 5)
        hold off
        plot(stim, 'y')
        hold on
        plot(i, stim(i), 'or', 'MarkerFaceColor', 'r')
    end
end



%% Guillaume's data:

addpath('/home/ljp/Science/Guillaume/Thermotaxis/Datasets')
h5disp('20171115_Run06_rp_Tset=27.h5')



clear; close all; clc



%% Defining focus:

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
dffanim = 0;
if dffanim
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:3000
        pause(0.000001)
        plot(dff(:, i))
        title(['Iteration ', num2str(i)])
        axis([1, 3000, -0.2, 1])
    end
end



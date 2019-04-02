clear; close all; clc



%% Import HDF5 file

cd '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/2018-05-24/Run 08/Analysis/HDF5'
hippo = 'h5hippo.h5';
% stim = h5read(hippo, '/Data/Stimulus');
hippotemp = H5F.open(hippo);
name = H5F.get_name(hippotemp);



%% Trying function h5switch:

nfolder = struct;
nfolder.path = '/home/ljp/Science/Hippolyte';
nfolder.name = 'newHDF5';
nfolder.filename = 'testfile1.h5';

information = struct;
% information.description = 'test to see.';
information.autodescription = true;
information.date = '21/10/2017';
information.run = 6;
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
overwrite.file = true;

h5switch(hippo, nfolder, information, overwrite);



%% Opening new HDF5 file:

hdf5test = fullfile(nfolder.path, nfolder.name, nfolder.filename);
h5disp(hdf5test)

dff = h5read(hdf5test, '/Data/Brain/Analysis/DFF');




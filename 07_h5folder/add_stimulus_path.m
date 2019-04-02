clear; close all; clc


% This program adds the path to the stimulus in the metadata for all HDF5 files
dataset_path = '/home/ljp/Science/Hippolyte/ALL_DATASETS';
dirdata = dir(dataset_path);

for i = 1:length(dirdata)
    ntemp = dirdata(i).name;
    npath = fullfile(dataset_path, ntemp);
    if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
        h5writeatt(npath, '/Metadata', 'Stimulus path', '/Data/Stimulus/vestibular1/motorAngle')
    elseif regexp(ntemp, '201\d{5}_Run\d{2}_rp_Tset=\d{1,}.h5')
        h5writeatt(npath, '/Metadata', 'Stimulus path', '/Data/Stimulus/RandomPulses/Trace')
    elseif regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}_a.h5')
        h5writeatt(npath, '/Metadata', 'Stimulus path', '/Data/Stimulus/auditory1/acousticPulse')
    end
end
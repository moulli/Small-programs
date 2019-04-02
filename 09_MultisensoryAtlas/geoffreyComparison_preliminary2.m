clear; close all; clc
addpath(genpath('/home/ljp/Programs'));
addpath(genpath('/home/ljp/Science/Hippolyte/Small-programs'))



%% Loading 3 step-stimulus files from Geoffrey & Guillaume:

% Geoffrey:
h1 = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run08';
h2 = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run11';
h3 = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles/2018-05-24Run15';
% Guillaume:
h4 = '/home/ljp/Science/Guillaume/Thermotaxis/Datasets/20171115_Run06_rp_Tset=27.h5';
h5 = '/home/ljp/Science/Guillaume/Thermotaxis/Datasets/20171116_Run05_rp_Tset=27.h5';
h6 = '/home/ljp/Science/Guillaume/Thermotaxis/Datasets/20180110_Run03_rp_Tset=28.h5';
% Building H:
H = {h1, h2, h3, h4, h5, h6};
nh = length(H);



%% Building grid from HDF5s:

method = 'Correlation analysis';
increment = 0.01;
zgrid = ZBraingrid(method, increment);
for i = [1:nh, 3]
    
    % Creating structure:
    attributes_in = struct;
    attributes_in.path = H{i};
    find_201 = strfind(H{i}, '201');
    attributes_in.name = H{i}(find_201(end):end);
    if i <= 3
        attributes_in.comment = 'vestibular + step';
    else
        attributes_in.comment = 'thermotaxis + random pulses';
    end
    
    % Loading data:
    attributes_in.coordinates = h5read(H{i}, '/Data/Brain/ZBrainCoordinates');
    dfftemp = h5read(H{i}, '/Data/Brain/Analysis/DFF');
    if i <= 3
        stimtemp = h5read(H{i}, '/Data/Stimulus/vestibular1/motorAngle');
    else
        stimtemp = h5read(H{i}, '/Data/Stimulus/RandomPulses/Trace');
    end
    stimtemp = reshape(stimtemp, 1, length(stimtemp));
    % Correlation coefficient:
    attributes_in.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
    
    
    % Adding example to zbraingrid object:
    addDataset(zgrid, attributes_in);
    disp(zgrid)
    
end
cleanDuplicates(zgrid)



%% Plotting:

ridval = 0.1;
plotAll(zgrid, 'ridvalues', [-ridval, ridval], 'MarkerSize', 40)
plotAll(zgrid, 'ridvalues', [-ridval, ridval], 'MarkerSize', 40, 'intercept', true)
plotSome(zgrid, 'vestibular', 'ridvalues', [-ridval, ridval], 'MarkerSize', 40)
plotSome(zgrid, 'vestibular', 'ridvalues', [-ridval, ridval], 'MarkerSize', 40, 'intercept', true)
%plotSome(zgridV, 'thermotaxis', 'ridvalues', [-ridval, ridval], 'MarkerSize', 40)







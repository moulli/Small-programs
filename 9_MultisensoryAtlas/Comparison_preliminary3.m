clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/Small-programs'))



%% Building structure:

method = 'Correlation analysis';
increment = 0.005;
zgrid = ZBraingrid(method, increment);



%% Scrolling through Geoffrey's and Guillaume's file, to add to zgrid:

geoffrey = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles';
dirgeo = dir(geoffrey);
for i = 1:length(dirgeo)
    ntemp = dirgeo(i).name;
    ptemp = fullfile(geoffrey, ntemp);
    if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
        % Info on dataset:
        stemp = struct;
        stemp.name = ntemp;
        stemp.path = ptemp;
        stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 sensory type')) + " + " + ...
                               string(h5readatt(ptemp, '/Metadata', 'Stimulus --> vestibular1 stimulus type'));
        % Data from dataset:
        stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
        dfftemp = h5read(ptemp, '/Data/Brain/Analysis/DFF');
        stimtemp = h5read(ptemp, '/Data/Stimulus/vestibular1/motorAngle');
        stemp.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
        addDataset(zgrid, stemp);
        disp(zgrid)
    end
end
        
guillaume = '/home/ljp/Science/Guillaume/Thermotaxis/Datasets';
dirgui = dir(guillaume);
for i = 1:length(dirgui)
    ntemp = dirgui(i).name;
    ptemp = fullfile(guillaume, ntemp);
    if regexp(ntemp, '201\d{5}_Run\d{2}_rp_Tset=\d{1,}.h5')
        % Info on dataset:
        stemp = struct;
        stemp.name = ntemp;
        stemp.path = ptemp;
        temperature = str2double(ntemp(24:end-3));
        if temperature <= 20
            templevel = "cold";
        elseif temperature >= 30
            templevel = "hot";
        else
            templevel = "neutral";
        end
        stemp.comment = string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses sensory type')) + " + " + ...
                               string(h5readatt(ptemp, '/Metadata', 'Stimulus --> RandomPulses stimulus type') + " + " + ...
                               string(temperature) + " degrees + " + templevel);
        % Data from dataset:
        stemp.coordinates = h5read(ptemp, '/Data/Brain/ZBrainCoordinates');
        dfftemp = h5read(ptemp, '/Data/Brain/Analysis/DFF');
        stimtemp = h5read(ptemp, '/Data/Stimulus/RandomPulses/Trace');
        stimtemp = stimtemp';
        stemp.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
        addDataset(zgrid, stemp);
        disp(zgrid)
    end
end



% %% Cleaning duplicates if necessary (should not be):
% 
% cleanDuplicates(zgrid);










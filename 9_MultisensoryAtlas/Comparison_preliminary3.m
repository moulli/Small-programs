clear; close all; clc
addpath(genpath('/home/ljp/Programs'))
addpath(genpath('/home/ljp/Science/Hippolyte/Small-programs'))



%% Building structure:

method = 'Correlation analysis';
increment = 0.005;
zgridGeo005 = ZBraingrid(method, increment);
zgridGui005 = ZBraingrid(method, increment);
 


%% Scrolling through Geoffrey's file, to add to zgrid:

geoffrey = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data/HDF5normalizedFiles';
dirgeo = dir(geoffrey);
for i = 1:length(dirgeo)
    ntemp = dirgeo(i).name;
    ptemp = fullfile(geoffrey, ntemp);
    if regexp(ntemp, '201\d(-\d{2}){2}Run\d{2}.h5')
        try
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
            % Normalizing:
            dfftemp = (dfftemp - mean(dfftemp, 2)) ./ std(dfftemp, [], 2);
            stimtemp = (stimtemp - mean(stimtemp)) ./ std(stimtemp);
            % Computing correlation and adding file:
            stemp.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
            addDataset(zgridGeo005, stemp);
            disp(zgridGeo005)
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
end



%% Cleaning duplicates if necessary (should not be):

cleanDuplicates(zgridGeo005);



%% Saving file:

pathcreatedGeo005 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo005.mat');
save(pathcreatedGeo005, 'zgridGeo005')



%% Saving with increments 0.01 and 0.05:

zgridGeo01 = downIncrement(zgridGeo005, 0.01);
pathcreatedGeo01 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo01.mat');
save(pathcreatedGeo01, 'zgridGeo01')

zgridGeo05 = downIncrement(zgridGeo005, 0.05);
pathcreatedGeo05 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGeo05.mat');
save(pathcreatedGeo05, 'zgridGeo05')



%% Scrolling through Guillaume's file, to add to zgrid:

guillaume = '/home/ljp/Science/Guillaume/Thermotaxis/Datasets';
dirgui = dir(guillaume);
for i = 1:length(dirgui)
    ntemp = dirgui(i).name;
    ptemp = fullfile(guillaume, ntemp);
    if regexp(ntemp, '201\d{5}_Run\d{2}_rp_Tset=\d{1,}.h5')
        try
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
            stimtemp = reshape(stimtemp, 1, length(stimtemp));
            % Normalizing:
            dfftemp = (dfftemp - mean(dfftemp, 2)) ./ std(dfftemp, [], 2);
            stimtemp = (stimtemp - mean(stimtemp)) ./ std(stimtemp);
            % Computing correlation and adding file:
            stemp.correlation = sum(dfftemp .* stimtemp, 2) ./ sqrt(sum(dfftemp.^2, 2) .* sum(stimtemp.^2, 2));
            addDataset(zgridGui005, stemp);
            disp(zgridGui005)
        catch
            fprintf('Problem with HDF5, moving on to the next one. \n');
        end
    end
end



%% Cleaning duplicates if necessary (should not be):

cleanDuplicates(zgridGui005);



%% Saving file:

pathcreatedGui005 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGui005.mat');
save(pathcreatedGui005, 'zgridGui005')



%% Saving with increments 0.01 and 0.05:

zgridGui01 = downIncrement(zgridGui005, 0.01);
pathcreatedGui01 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGui01.mat');
save(pathcreatedGui01, 'zgridGui01')

zgridGui05 = downIncrement(zgridGui005, 0.05);
pathcreatedGui05 = fullfile('/home/ljp/Science/Hippolyte', 'zgridGui05.mat');
save(pathcreatedGui05, 'zgridGui05')



%% Assiging new datasets:

zgrid = zgridGeo005 + zgridGui005;



%% Studying correlation distibutions:

corre = cell(5, 1);
segments = ["step", "sine", "hot", "cold", "neutral"];
for i = 1:5
    zgridtemp = choseSubset(zgrid, segments(i));
    lzgrid = length(zgridtemp);
    corretemp = zeros(lzgrid, 2);
    for j = 1:lzgrid
        corretemp(j, 1) = mean(zgridtemp.Zcorvect{j});
        corretemp(j, 2) = var(zgridtemp.Zcorvect{j});
    end
    corre{i} = corretemp;
end

cormat = zeros(5, 4);
for i = 1:5
    cormat(i, 1) = mean(corre{i}(:, 1));
    cormat(i, 2) = var(corre{i}(:, 1)); 
    cormat(i, 3) = mean(corre{i}(:, 2));
    cormat(i, 4) = var(corre{i}(:, 2)); 
end
titles = ["Mean correlation across datasets", "Mean correlation variance across datasets", ...
          "Correlation variance across datasets", "Correlation variance variance across datasets"];
subsets = ["Vestibular step", "Vestibular sine", "Thermotaxis hot", "Thermotaxis cold", "Thermotaxis neutral"];
figure
for i = 1:4
    subplot(4, 1, i)
    cortemp = cormat(:, i);
    plot(cortemp, ':o', 'MarkerFaceColor', [0, 0, 0])
    title(titles(i), 'Interpreter', 'latex')
    difftemp = max(cortemp) - min(cortemp);
    axis([0.9, 5.1, min(cortemp)-0.1*difftemp, max(cortemp)+0.1*difftemp])
    xticks(1:5)
    xticklabels(subsets)
    grid on
end










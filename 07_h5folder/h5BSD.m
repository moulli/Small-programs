clear; close all; clc



%% Including a few tools:

addpath(genpath('/home/ljp/Programs/'));
addpath(genpath('/home/ljp/Science/Hippolyte/Small-programs'));



%% Indicating old folder path, new folder path, and new folder name:

folder_path = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data';



%% Creating loop for each file in folder:

% Getting file with already created HDF5, or building one:
pathcreated = fullfile(folder_path, 'HDF5createdBSD.mat');
try
    load(pathcreated)
catch
    HDF5createdBSD = zeros(0, 2);
    save(pathcreated, 'HDF5createdBSD')
end

% Creating structure with files:
all_files = dir(folder_path);

% Launching loop:
for i = 1:length(all_files)
    
    filetemp = all_files(i).name;
    % If file is 10 characters and starts with 201, checking runs:
    if length(filetemp) == 10 && startsWith(filetemp, '201')
        
        % Storing data:
        root = folder_path;
        if endsWith(folder_path, 'Data')
            root = root(1:end-4);
        end
        study = '';
        date = filetemp;
        
        % Checking runs:
        all_runs = dir([folder_path, '/', filetemp]);
        for j = 1:length(all_runs)
            runtemp = all_runs(j).name;
            % If date and run already processed, going on:
            if ~isempty(HDF5createdBSD)
                if sum(filetemp == HDF5createdBSD(:, 1)) > 0 && sum(runtemp == HDF5createdBSD(:, 2)) > 0
                    fprintf('HDF5 for date %s and %s already created, moving on to next one. \n', filetemp, runtemp);
                    continue
                end
            end
            % If file starts with 'Run', creating focus:
            if startsWith(runtemp, 'Run')
                
                % Checking if 'Analysis' folder exists:
                if ~exist([folder_path, '/', filetemp, '/', runtemp, '/', 'Analysis'], 'dir')
                    continue
                end
                
                % Focus:
                run = runtemp;
                try
                    warning('off')
                    F = NT.Focus(root, study, date, run);
                    warning('on')
                catch
                    continue
                end
                
                
                % Adding data to already existing HDF5:
                
                    % Checking if HDF5 already exists (necessary):
                    try
                        lstemp = ls(F.dir('HDF5'));
                    catch
                        fprintf('No HDF5 file in the focus, moving on to next one. \n');
                        continue
                    end
                    hdf5path = fullfile(F.dir('HDF5'), [filetemp, '(', erase(runtemp, ' '), ').h5']);
                    dfftemp = cast(h5read(hdf5path, '/Data/Values'), 'double');
                    dfftemp(abs(dfftemp) < 0.000001) = 0.000001; 
                    times = cast(h5read(hdf5path, '/Data/Times'), 'double');
                    nneu = size(dfftemp, 1);
                    
                
                    % Informing on which file we are doing:
                    fprintf('\nCreating BSD and Regression label in HDF5 for %s, date: %s. \n', runtemp, filetemp);
                    
                    % BSD program:
                        % Defining structure:
                        Oalg = struct; 
                        Oalg.Time = size(dfftemp, 2); 
                        Oalg.dt = mean(gradient(times)); 
                        Oalg.nNeurons = nneu; 
                        Oalg.adaptive = 0; 
                        % Running algorithm:
                        [Nspikes, ~, ~, Pphys, ~] = pBSD(dfftemp, Oalg);
                        % Adding data to HDF5:
                        h5create(hdf5path, '/Data/Brain/Analysis/InferredSpikes', size(Nspikes), 'Datatype', 'single');
                        h5write(hdf5path, '/Data/Brain/Analysis/InferredSpikes', cast(Nspikes, 'single'));
                        h5create(hdf5path, '/Data/Brain/Analysis/ThresholdedSpikes', size(Nspikes), 'Datatype', 'uint8');
                        h5write(hdf5path, '/Data/Brain/Analysis/ThresholdedSpikes', cast(1 * (Nspikes > Pphys.threshold'), 'uint8'));
                    
                    
                    % Regression:
                        % Launching h5stimreg:
                        Params.period = 1;
                        [out, ~, ~] = h5stimreg(hdf5path, Params);
                        % Quantile and parameters:
                        q = 10;
                        coef = sqrt(sum((out.coef.^2), 2));
                        fstat = out.Fstat;                    
                        % Coefficients:
                        [maxcoef, indcoef] = sort(coef, 'descend');
                        indcK = sort(indcoef(1:ceil(nneu/q)));
                        % F-stat:
                        [maxfstat, indfstat] = sort(fstat, 'descend');
                        indfK = sort(indfstat(1:ceil(nneu/q)));
                        % Taking indices that are in both vectors:
                        compa = sum(indcK == indfK', 2);
                        indfin = indcK(compa == 1);
                        % Building labels vector:
                        RegLabels = zeros(nneu, 1);
                        RegLabels(indfin) = 1;
                        % Adding data to HDF5:
                        h5create(hdf5path, '/Data/Brain/Analysis/RegressionLabels', size(RegLabels), 'Datatype', 'single');
                        h5write(hdf5path, '/Data/Brain/Analysis/RegressionLabels', cast(RegLabels, 'single'))
                        
                        
                    % Saving date and run in the HDF5created matrix:
                    HDF5createdBSD = [HDF5createdBSD; string(filetemp), string(runtemp)];
                    save(pathcreated, 'HDF5createdBSD')
                    
                    
               
            end
        end
        
    end
    
end
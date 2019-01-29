clear; close all; clc



%% Including a few tools:

addpath(genpath('/home/ljp/Programs/'));
addpath(genpath('/home/ljp/Science/Hippolyte/Small-programs'));



%% Indicating old folder path, new folder path, and new folder name:

old_folder_path = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data';
new_folder_path = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault';
new_folder_name =  'HDF5createdFiles';



%% Creating loop for each file in folder:

% Getting file with already created HDF5, or building one:
pathcreated = fullfile(new_folder_path, 'HDF5created.mat');
try
    load(pathcreated)
catch
    HDF5created = zeros(0, 2);
    save(pathcreated, 'HDF5created')
end

% Creating structure with files:
all_files = dir(old_folder_path);

% Launching loop:
for i = 1:length(all_files)
    
    filetemp = all_files(i).name;
    % If file is 10 characters and starts with 201, checking runs:
    if length(filetemp) == 10 && startsWith(filetemp, '201')
        
        % Storing data:
        root = old_folder_path;
        if endsWith(old_folder_path, 'Data')
            root = root(1:end-4);
        end
        study = '';
        date = filetemp;
        
        % Checking runs:
        all_runs = dir([old_folder_path, '/', filetemp]);
        for j = 1:length(all_runs)
            runtemp = all_runs(j).name;
            % If date and run already processed, going on:
            if ~isempty(HDF5created)
                if sum(filetemp == HDF5created(:, 1)) > 0 && sum(runtemp == HDF5created(:, 2)) > 0
                    fprintf('HDF5 for date %s and %s already created, moving on to next one. \n', filetemp, runtemp);
                    continue
                end
            end
            % If file starts with 'Run', creating focus:
            if startsWith(runtemp, 'Run')
                
                % Checking if 'Analysis' folder exists:
                if ~exist([old_folder_path, '/', filetemp, '/', runtemp, '/', 'Analysis'], 'dir')
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
                
                
                % Creating structures:
                
                    % Checking if HDF5 already exists (necessary):
                    try
                        lstemp = ls(F.dir('HDF5'));
                    catch
                        fprintf('No HDF5 file in the focus, moving on to next one. \n');
                        continue
                    end
                    
                
                    % Informing on which file we are doing:
                    fprintf('\nCreating new HDF5 for %s, date: %s. \n', runtemp, filetemp);

                    % Creating nfolder structure:
                    nfolder = struct;
                    nfolder.path = new_folder_path;
                    nfolder.name = new_folder_name;
                    nfolder.filename = [filetemp, erase(runtemp, ' ')];
                    % Deleting HDF5 if not finished:
                    if isfile(fullfile(nfolder.path, nfolder.name, nfolder.filename))
                        delete(fullfile(nfolder.path, nfolder.name, nfolder.filename))
                    end

                    % Creating overwrite structure:
                    overwrite = struct;
                    overwrite.folder = false;
                    overwrite.file = false;

                    % Creating information structure:
                    information = h5writeInfoStruct(F);
                    
                    
                % h5switch function:
                
                    h5switch(F, nfolder, information, overwrite);
                    

                % Saving date and run in the HDF5created matrix:
                
                    HDF5created = [HDF5created; string(filetemp), string(runtemp)];
                    save(pathcreated, 'HDF5created')
                    
               
            end
        end
        
    end
    
end









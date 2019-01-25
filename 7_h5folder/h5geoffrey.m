clear; close all; clc



%% Including a few tools:

addpath(genpath('/home/ljp/Programs/'));
addpath(genpath('/home/ljp/Science/Hippolyte/Small-programs'));



%% Indicating old folder path, new folder path, and new folder name:

old_folder_path = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault/Data';
new_folder_path = '/home/ljp/Science/GeoffreysComputer/Paper_Data/2018_Migault';
new_folder_name = 'DataNormalizedHDF5';



%% Creating loop for each file in folder:

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
                
                    % Informing on which file we are doing:
                    fprintf('\nCreating new HDF5 for %s, date: %s. \n', runtemp, filetemp);

                    % Creating nfolder structure:
                    nfolder = struct;
                    nfolder.path = new_folder_path;
                    nfolder.name = new_folder_name;
                    nfolder.filename = [filetemp, erase(runtemp, ' ')];

                    % Creating overwrite structure:
                    overwrite = struct;
                    overwrite.folder = false;
                    overwrite.file = false;

                    % Creating information structure:
                    information = h5writeInfoStruct(F);
%                     
%                     
%                 % h5switch function:
%                 
%                     h5switch(F, nfolder, information, overwrite);

%% ADD WHAT HAS ALREADY BEEN DONE!!!!
                    
               
            end
        end
        
    end
    
end









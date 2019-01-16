function h5switch(h5old_Geoffrey, nfolder, information, overwrite)



    %% Checking structures:
    
    % Indication:
    tic
    fprintf('\n\nFunction h5switch started\n\n');
    
    % Creation of silent boolean:
    if ~isfield(overwrite, 'folder'); overwrite.folder = false; end
    if ~isfield(overwrite, 'file'); overwrite.file = false; end
    
    % Structure: folder:
    if ~isfield(nfolder, 'path'); error('Please provide folder path.'); end
    if ~isfield(nfolder, 'name'); error('Please provide folder name.'); end
    if ~isfield(nfolder, 'filename'); error('Please provide file name.'); end
    
    % Structure: information:
    if ~isfield(information, 'description'); information.description = ''; end
    if ~isfield(information, 'autodescription'); information.description = false; end
    if ~isfield(information, 'date'); error('Please provide date.'); end
    if ~isfield(information, 'run'); error('Please provide run.'); end
    if ~isfield(information, 'line'); error('Please provide line.'); end
    if ~isfield(information, 'age'); error('Please provide age.'); end
    if ~isfield(information, 'autoID'); information.autoID = false; end
    if ~isfield(information, 'ID') && ~isfield(information, 'autoID')
        error('Please provide ID.')
    elseif ~isfield(information, 'ID') && isfield(information, 'autoID')
        information.ID = erase([information.date, num2str(information.run)], [" ", "/"]);
    end
    if isfield(information, 'stimulus')
        if ~isfield(information.stimulus, 'name'); error('Please provide name of stimulus/stimuli.'); else; n1 = size(information.stimulus.name); end
        if ~isfield(information.stimulus, 'sensorytype'); error('Please provide sensory type of stimulus/stimuli.'); else; n2 = size(information.stimulus.sensorytype); end
        if ~isfield(information.stimulus, 'stimtype'); error('Please provide type/types of stimulation.'); else; n3 = size(information.stimulus.stimtype); end
        if ~isfield(information.stimulus, 'frequency'); error('Please provide frequency of stimulus/stimuli.'); else; n4 = size(information.stimulus.frequency); end
        if ~isequal([n1; n2; n3; n4], n1 .* ones(4, 2)); error('Please provide stimulus information with same sizes.'); end
    end
    if isfield(information, 'behaviour')
        if ~isfield(information.behaviour, 'name'); error('Please provide name of behaviour/s.'); else; nb1 = size(information.behaviour.name); end
        if ~isfield(information.behaviour, 'type'); error('Please provide type of behaviour/s.'); else; nb2 = size(information.behaviour.type); end
        if ~isfield(information.behaviour, 'frequency'); error('Please provide frequency of behaviour/s.'); else; nb3 = size(information.behaviour.frequency); end
        if ~isequal([nb1; nb2; nb3], nb1 .* ones(3, 2)); error('Please provide stimulus information with same sizes.'); end
    end


    
    %% Creating folder for new HDF5 file:
    
    % Changing directory:
    curdir = pwd;
    cd(nfolder.path)
    folderall = fullfile(nfolder.path, nfolder.name);
    % Folder:
    [~, ~, msgID] = mkdir(nfolder.name);
    switch msgID
        case 'MATLAB:MKDIR:DirectoryExists'
            fprintf('Folder name already exists. \n');
            if overwrite.folder
                prompt = 'Are you sure you want to overwrite whole folder? Y/N [Y]: ';
                inStr = input(prompt,'s');
                switch lower(inStr)
                    case 'y'
                        fprintf('Overwriting folder %s. \n', folderall);
                        rmdir(nfolder.name, 's')
                        mkdir(nfolder.name)
                        fprintf('Folder deleted and new folder created at %s. \n', folderall);
                    otherwise
                        error('Please launch function again with appropriate parameters.')
                end
            else
                fprintf('Keeping already existing folder. \n');
            end
        otherwise
            mkdir(nfolder.name)
            fprintf('Folder created at %s. \n', folderall);
    end
    cd(folderall)
    
    
    
    %% Creating new HDF5 file:
    
    % Checking if file already exists:
    switch isfile(nfolder.filename)
        case 1
            fprintf('File name already exists in folder. \n');
            if overwrite.file
                prompt = 'Are you sure you want to overwrite file? Y/N [Y]: ';
                inStr = input(prompt,'s');
                switch lower(inStr)
                    case 'y'
                        fprintf('Overwriting file %s. \n', nfolder.filename);
                        delete(nfolder.filename)
                        fprintf('File %s deleted and new file created. \n', nfolder.filename);
                    otherwise
                        error('Please launch function again with appropriate parameters.')
                end 
            else
                error('Conflict between names.')
            end
        otherwise
            fprintf('File %s created. \n', nfolder.name);
    end
    
    
    
    %% Metadata:
    
    % Creation:
    h5create(nfolder.filename, '/Metadata', [1, 1], 'Datatype', 'double');
    h5write(nfolder.filename, '/Metadata', 1)
    % Experiment:
    h5writeatt(nfolder.filename, '/Metadata', 'Experiment date (jj/mm/dddd)', information.date)
    h5writeatt(nfolder.filename, '/Metadata', 'Experiment run', uint8(information.run))
    % Fish:
    h5writeatt(nfolder.filename, '/Metadata', 'Fish line', information.line)
    h5writeatt(nfolder.filename, '/Metadata', 'Fish age (dpf)', uint8(information.age))
    h5writeatt(nfolder.filename, '/Metadata', 'Fish ID', information.ID)
    % Stimulus:
    for i = 1:max(n1)
        h5writeatt(nfolder.filename, '/Metadata', ['Stimulus --> ', information.stimulus.name{i}, ' sensory type'], information.stimulus.sensorytype{i})
        h5writeatt(nfolder.filename, '/Metadata', ['Stimulus --> ', information.stimulus.name{i}, ' stimulus type'], information.stimulus.stimtype{i})
        h5writeatt(nfolder.filename, '/Metadata', ['Stimulus --> ', information.stimulus.name{i}, ' frequency (Hz)'], information.stimulus.frequency{i})
    end
    % Behaviour:
    for i = 1:max(nb1)
        h5writeatt(nfolder.filename, '/Metadata', ['Behaviour --> ', information.behaviour.name{i}, ' type'], information.behaviour.type{i})
        h5writeatt(nfolder.filename, '/Metadata', ['Behaviour --> ', information.behaviour.name{i}, ' frequency (Hz)'], information.behaviour.frequency{i})
    end
    % File:
    if isfield(information, 'creation'); h5writeatt(nfolder.filename, '/Metadata', 'File creation', information.creation); end
    if isfield(information, 'programname'); h5writeatt(nfolder.filename, '/Metadata', 'File program name', information.programname); end
    if isfield(information, 'programhash'); h5writeatt(nfolder.filename, '/Metadata', 'File program hash', information.programhash); end
    
    
    
    
    %% Data:
    
    % Getting path to HDF5 file:
    h5old_G = [curdir, '/', h5old_Geoffrey];
    % Times:
    Times = h5read(h5old_G, '/Data/Times');
    h5create(nfolder.filename, '/Data/Brain/Time', size(Times), 'Datatype', 'double');
    h5write(nfolder.filename, '/Data/Brain/Time', Times)
    h5writeatt(nfolder.filename, '/Data/Brain/Time', 'unit', 's')    
    % Time delays:
    TimeDelays = h5read(h5old_G, '/Data/TimeDelays');
    h5create(nfolder.filename, '/Data/Brain/TimeDelays', size(TimeDelays), 'Datatype', 'double');
    h5write(nfolder.filename, '/Data/Brain/TimeDelays', TimeDelays)
    h5writeatt(nfolder.filename, '/Data/Brain/TimeDelays', 'unit', 's') 
    % Coordinates:
    Coordinates = h5read(h5old_G, '/Data/Coordinates');
    h5create(nfolder.filename, '/Data/Brain/Coordinates', size(Coordinates), 'Datatype', 'double');
    h5write(nfolder.filename, '/Data/Brain/Coordinates', Coordinates)
    h5writeatt(nfolder.filename, '/Data/Brain/Coordinates', 'unit', 'mm') 
    h5writeatt(nfolder.filename, '/Data/Brain/Coordinates', 'space', 'RAS') 
    % Reference coordinates:
    RefCoordinates = h5read(h5old_G, '/Data/RefCoordinates');
    h5create(nfolder.filename, '/Data/Brain/RefCoordinates', size(RefCoordinates), 'Datatype', 'double');
    h5write(nfolder.filename, '/Data/Brain/RefCoordinates', RefCoordinates)
    h5writeatt(nfolder.filename, '/Data/Brain/RefCoordinates', 'unit', 'mm') 
    h5writeatt(nfolder.filename, '/Data/Brain/RefCoordinates', 'space', 'RAS') 
    h5writeatt(nfolder.filename, '/Data/Brain/RefCoordinates', 'reference brain', 'zbrain atlas')  
    % Labels:
    Labels = h5read(h5old_G, '/Data/Labels');
    h5create(nfolder.filename, '/Data/Brain/Labels', size(Labels), 'Datatype', 'double');
    h5write(nfolder.filename, '/Data/Brain/Labels', Labels)
    h5writeatt(nfolder.filename, '/Data/Brain/Labels', 'origin', 'zbain atlas') 
    % DFF:
    Values = h5read(h5old_G, '/Data/Values');
    h5create(nfolder.filename, '/Data/Brain/Analysis/Values', size(Values), 'Datatype', 'double');
    h5write(nfolder.filename, '/Data/Brain/Analysis/Values', Values)
    % Stimulus:
    Stimulus = h5read(h5old_G, '/Data/Stimulus');
    h5create(nfolder.filename, ['/Data/Stimulus/', information.stimulus.name{i}, '/motorAngle'], size(Stimulus), 'Datatype', 'double')
    h5write(nfolder.filename, ['/Data/Stimulus/', information.stimulus.name{i}, '/motorAngle'], Stimulus)
    h5writeatt(nfolder.filename, ['/Data/Stimulus/', information.stimulus.name{i}, '/motorAngle'], 'unit', 'degrees')
    
    
    
    %% Description:
    
    h5create(nfolder.filename, '/Description', [1, 1], 'Datatype', 'double');
    if ~isempty(information.description)
        h5write(nfolder.filename, '/Description', 1)
        h5writeatt(nfolder.filename, '/Description', 'Description', information.description)
    elseif isempty(information.description) && ~information.autodescription
        h5write(nfolder.filename, '/Description', 0)
    else 
        autodes = [information.date, ', Run ', num2str(information.run)];
        h5writeatt(nfolder.filename, '/Description', 'Description', autodes)
    end
    
    
    
    %% Changing back to initial directory:
    
    cd(curdir)
    fprintf('\nFunction h5switch ended in %.3f seconds. \n\n\n', toc);
    


end
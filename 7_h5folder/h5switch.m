function h5switch(F, nfolder, information, overwrite)



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
    if ~isfield(information, 'rate'); error('Please provide acquisition rate.'); end
    if ~isfield(information, 'layers'); error('Please provide number of layers.'); end
    if ~isfield(information, 'increment'); error('Please provide increment distance.'); end
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
                inStr = input(prompt, 's');
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
                inStr = input(prompt, 's');
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
    h5writeatt(nfolder.filename, '/Metadata', 'Experiment date (dddd/mm/jj)', information.date)
    h5writeatt(nfolder.filename, '/Metadata', 'Experiment run', uint8(information.run))
    % Acquisition:
    h5writeatt(nfolder.filename, '/Metadata', 'Acquisition rate (Hz)', information.rate)
    h5writeatt(nfolder.filename, '/Metadata', 'Acquisition number of layers', uint8(information.layers))
    h5writeatt(nfolder.filename, '/Metadata', 'Acquisition interlayer distance', information.increment)
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
    
    % Getting path to HDF5 file from focus:
    filetemp = strsplit(ls(F.dir('HDF5')));
    for i = 1:length(filetemp)
        if isequal(filetemp{i}(end-2:end), '.h5')
            filetemp = filetemp{i};
            break
        end
    end
    h5old_G = [F.dir('HDF5'), '/', filetemp];
    % Times:
    h5fill(h5old_G, '/Data/Times', nfolder.filename, '/Data/Brain/Times', 'single', {'unit', 's'})
    % Time delays:
    h5fill(h5old_G, '/Data/TimeDelays', nfolder.filename, '/Data/Brain/TimeDelays', 'single', {'unit', 's'})
    % Coordinates:
    h5fill(h5old_G, '/Data/Coordinates', nfolder.filename, '/Data/Brain/Coordinates', 'single', {'unit', 'mm'; 'space', 'RAS'})
    % Reference coordinates:
    h5fill(h5old_G, '/Data/RefCoordinates', nfolder.filename, '/Data/Brain/RefCoordinates', 'single', {'unit', 'mm'; 'space', 'RAS'; 'reference brain', 'zbrain atlas'})
    % Labels:
    h5fill(h5old_G, '/Data/Labels', nfolder.filename, '/Data/Brain/Labels', 'single', {'origin', 'zbrain atlas'})
    % DFF:
    h5fill(h5old_G, '/Data/Values', nfolder.filename, '/Data/Brain/Analysis/DFF', 'single')
    % Stimulus:
    h5fill(h5old_G, '/Data/Stimulus', nfolder.filename, ['/Data/Stimulus/', information.stimulus.name{i}, '/motorAngle'], 'single', {'unit', 'degrees'})
    
    
    
    %% Rebuilding raw signal if necessary:
    
    % Indication:
    fprintf('\nMetadata and data filled in %.3f s; building raw data. \n', toc);
    % Create corrected stack:
    try 
        existcorrected = 1;
        m = Focused.Mmap(F, 'corrected');
        fprintf('Corrected stack already exists. \n');
    catch
        existcorrected = 0;
        fprintf('Building corrected stack. \n');
        m = Focused.Mmap(F, 'graystack');
        F.Analysis.Layers = m.Z;
        driftApply(F);
        m = Focused.Mmap(F, 'corrected');
        fprintf('Corrected stack built (%.3f s). \n', toc);
    end
    % Define raw signal matrix:
    sizeHDF5 = size(h5read(h5old_G, '/Data/TimeDelays'), 1);
    rawtrace = zeros(sizeHDF5, m.t);
    % Load neuron shape
    fprintf('Launching average computation for all neurons segmented. \n');
    neuCount = 1;
    layer = 1;
    for iz = m.Z
        segPath = F.dir('Segmentation');
        inputSeg = fullfile(segPath, [num2str(iz, '%02d') '.mat']); % iz = layer
        load(inputSeg, 'neuronShape');
        ns = neuronShape;
        % Average on pixel:
        for i = 1:length(ns)
            rawtrace(neuCount, :) = squeeze(mean(m(ns{i}, iz, :), 1))';
            if mod(i, 1000) == 0
                fprintf('Iteration %.0f out of %.0f, for layer %.0f out of %.0f, in %.3f seconds. \n', [i, length(ns), layer, length(m.Z), toc]);
            end
            neuCount = neuCount + 1;
        end
        layer = layer + 1;
    end
    % Adding raw signal in HDF5:
    h5create(nfolder.filename, '/Data/Brain/RawSignal', size(rawtrace), 'Datatype', 'single');
    h5write(nfolder.filename, '/Data/Brain/RawSignal', cast(rawtrace, 'single'))
    % Deleting corrected stack:
    if existcorrected == 0
        while 1
            prompt = 'Do you want to keep corrected.stack? Y/N [Y]: ';
            inStr = input(prompt, 's');
            switch lower(inStr)
                case 'y'
                    break
                case 'n'
                    fclose('all');
                    rmdir(F.dir('corrected'), 's')
                    break
            end
        end
    end
        
    
    
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



%% Fill HDF5 function:
function h5fill(hdf5name, hdf5pathold, hdf5filename, hdf5pathnew, datatype, hdf5attributes)
% Function that automatically fills new HDF5 with old HDF5. 
% Inputs:
%    -- hdf5name: name of old HDF5.
%    -- hdf5pathold: path to data in old HDF5.
%    -- hdf5filename: name of new HDF5.
%    -- hdf5pathnew: path to data in new HDF5.
%    -- datatype: type of data for new HDF5.
%    -- hdf5attributes: {n x 2} cell of attributes for new HDF5.    
    
    % Getting info from old HDF5:
    temp = h5read(hdf5name, hdf5pathold);
    % Creating new path in new HDF5 and filling it:
    h5create(hdf5filename, hdf5pathnew, size(temp), 'Datatype', datatype);
    h5write(hdf5filename, hdf5pathnew, cast(temp, datatype))
    % Adding attributes:
    if nargin == 6
        numatt = size(hdf5attributes, 1);
        if numatt ~= 0
            for i = 1:numatt
                h5writeatt(hdf5filename, hdf5pathnew, hdf5attributes{i, 1}, hdf5attributes{i, 2})
            end
        end
    end

end
    
    
    
    
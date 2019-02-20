function addDataset(obj, dataset_in)

%% Function that adds information on new data in the ZBraingrid class.
%
%  New data is provided through a structure, including new dataset
%  coordinates and correlation coefficients. 
%
%
%% Inputs:
%
%  --obj: references the object this methods is attached to.
%  --dataset_in: structure containing fields 'name', 'path', 'comment'
%    (optional), 'coordinates', and 'correlation'. 'comment' should state
%    the type of stimulus provided (e.g. vestibular + step, or thermotaxis
%    + random pulses). 'coordinates' should have as many rows as there are
%    neurons in the dataset, and 3 columns for the 3 dimensions.
%    'correlation' should be a vector, with as many values as there are
%    neurons.



    %% Initialization:
    
    % Indication:
    tic
    if size(obj.Zneurons, 4) > 1 || sum(obj.Zneuron_number(:)) ~= 0 
        fprintf('\nLaunching function addDataset, attribute of ZBraingrid class. %.0f dataset(s) already added.\n', size(obj.Zneurons, 4));
    else
        fprintf('\nLaunching function addDataset, attribute of ZBraingrid class. First dataset.\n');
    end
    % Checking entering structure fields:
    if ~isfield(dataset_in, 'name'); error('Please provide dataset name.'); end
    if ~isfield(dataset_in, 'path'); error('Please provide dataset path.'); end
    if ~isfield(dataset_in, 'comment'); dataset_in.comment = 'None'; end
    if ~isfield(dataset_in, 'coordinates'); error('Please provide coordinates for each neuron.'); end
    if ~isfield(dataset_in, 'correlation'); error('Please provide correlation for each neuron.'); end
    % Checking dimensions:
    coord_in = dataset_in.coordinates;
    cor_in = dataset_in.correlation;
    [mcoord, ncoord] = size(coord_in);
    if ncoord ~= 3
        error('Coordinates matrix should have 3 columns for 3 dimensions.')
    elseif ~isvector(cor_in)
        error('Correlation should be a vector.')
    elseif mcoord ~= length(cor_in)
        error('Coordinates should have as many rows as there are correlations.')
    else 
        cor_in = reshape(cor_in, length(cor_in), 1);
    end
    % Adding data to properties:
    % Name:
    if isempty(obj.names)
        obj.names = {string(dataset_in.name)};
    else
        obj.names = [obj.names; string(dataset_in.name)];
    end
    % Path:
    if isempty(obj.paths)
        obj.paths = {string(dataset_in.path)};
    else
        obj.paths = [obj.paths; string(dataset_in.path)];
    end
    % Comment:
    if isempty(obj.comments)
        obj.comments = {string(dataset_in.comment)};
    else
        obj.comments = [obj.comments; string(dataset_in.comment)];
    end
    % Correlation vector:
    obj.Zcorvect = [obj.Zcorvect; cor_in];
    
    
    
    %% Filling new data to Zneurons, Zcorrelations & Zneuron_number:
    
    % Assigning new cell for Zneurons, new matrices for Zcorrelations & Zneuron_number, with new room:
    lgrid_in = size(obj.Zneurons);
    Zneurons_temp = cell(lgrid_in(1:3));
    Zcorrelations_temp = zeros(lgrid_in(1:3));
    Zneuron_number_temp = zeros(lgrid_in(1:3));
    % Algorithm:
    for ix = 1:lgrid_in(1)
        for iy = 1:lgrid_in(2)
            for iz = 1:lgrid_in(3)
                % Getting neurons in this part of grid:
                xtemp = (obj.xgrid(ix) <= coord_in(:, 1) & coord_in(:, 1) < obj.xgrid(ix+1));
                ytemp = (obj.ygrid(iy) <= coord_in(:, 2) & coord_in(:, 2) < obj.ygrid(iy+1));
                ztemp = (obj.zgrid(iz) <= coord_in(:, 3) & coord_in(:, 3) < obj.zgrid(iz+1));
                Ttemp = find(xtemp & ytemp & ztemp);
                % Now updating zbraingrid object:
                if isempty(Ttemp)
                    Zneurons_temp{ix, iy, iz} = [];
                else
                    Zneurons_temp{ix, iy, iz} = Ttemp;
                end
                cortemp = mean(cor_in(Ttemp));
                if isnan(cortemp); cortemp = 0; end
                Zcorrelations_temp(ix, iy, iz) = cortemp;
                Zneuron_number_temp(ix, iy, iz) = length(Ttemp);
            end
        end
        % Indication:
        if mod(ix, 4) == 0
            fprintf('For-loop %.2f %% completed in %.2f seconds.\n', [100*ix/lgrid_in(1), toc]);
        end
    end    
    % Adding new dataset to object:
    if size(obj.Zneurons, 4) > 1 || sum(obj.Zneuron_number(:)) ~= 0 
        obj.Zneurons = cat(4, obj.Zneurons, Zneurons_temp);
        obj.Zcorrelations = cat(4, obj.Zcorrelations, Zcorrelations_temp);
        obj.Zneuron_number = cat(4, obj.Zneuron_number, Zneuron_number_temp);
    else
        obj.Zneurons = Zneurons_temp;
        obj.Zcorrelations = Zcorrelations_temp;
        obj.Zneuron_number = Zneuron_number_temp;
    end
    obj.gridsize(4) = obj.gridsize(4) + 1;
    % Indication:
    fprintf('Function addDataset ended in %.2f seconds.\n', toc);
    
    
end
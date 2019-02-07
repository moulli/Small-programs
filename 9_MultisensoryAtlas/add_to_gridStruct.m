function gridStruct_out = add_to_gridStruct(gridStruct_in, name_in, zcoord_in, cor_in)

%% Function that adds information on new data in the gridStruct structure.
%
%  New data is provided through inputs cor_in and zcoord_in which represent
%  respectively correlation and coordinates. Then, grid is analyzed to
%  verify it includes all the new values (grid is extended otherwise), and
%  new information is computed.
%
%
%% Inputs:
%
%  --gridStruct_in: structure containing fields 'method', 'names', 'xgrid', 
%    'ygrid', 'zgrid', 'Tgrid' and 'Cgrid'. The first field is the method 
%    employed for the analysis, second field is the name of the HDF5 
%    already computed, third, fourth and fifth fields are information on 
%    the grid, the sixth one if a cell containing for each HDF5 all the 
%    neurons for each part of the grid, and seventh field is the mean
%    correlation for each part of the grid.
%  --name_in: string, name of the new file added.
%  --zcoord_in: (n x 3) matrix, n being the number of neurons, and
%    coordinates should be taken from the ZBrain Atlas.
%  --cor_in: array of length n, with correlation coefficient for each
%    neuron provided.
%
%
%% Output:
%
%  --gridStruct_out: modified gridStruct_in with information on new file.



    %% Initialization:
    
    % Checking entering structure fields:
    if ~isfield(gridStruct_in, 'method') || ~isfield(gridStruct_in, 'names') || ...
       ~isfield(gridStruct_in, 'xgrid') || ~isfield(gridStruct_in, 'ygrid') || ...
       ~isfield(gridStruct_in, 'zgrid') || ~isfield(gridStruct_in, 'Tgrid') || ...
       ~isfield(gridStruct_in, 'Cgrid')
        error('Please provide gridStruct_in with the required fields.')
    end
    % Checking entering structure interior dimensions:
    lnames_in = length(gridStruct_in.names);
    if ~isvector(gridStruct_in.xgrid); error('xgrid field must be a vector.'); end
    if ~isvector(gridStruct_in.ygrid); error('ygrid field must be a vector.'); end
    if ~isvector(gridStruct_in.zgrid); error('zgrid field must be a vector.'); end
    if length(size(gridStruct_in.Tgrid)) ~= 4
        error('Tgrid field must be a 4 dimension matrix.')
    end
    lTgrid_in = size(gridStruct_in.Tgrid, 4);
    if lnames_in ~= lTgrid_in; error('Number of names does not correspond to number of examples in Tgrid field.'); end
    if length(size(gridStruct_in.Cgrid)) ~= 4
        error('Cgrid field must be a 4 dimension matrix.')
    end
    lCgrid_in = size(gridStruct_in.Cgrid, 4);
    if lnames_in ~= lCgrid_in; error('Number of names does not correspond to number of examples in Cgrid field.'); end
    % Duplicating data from structure:
    method = gridStruct_in.method;
    names = gridStruct_in.names;
    xgrid = gridStruct_in.xgrid;
    ygrid = gridStruct_in.ygrid;
    zgrid = gridStruct_in.zgrid;
    Tgrid = gridStruct_in.Tgrid;
    Cgrid = gridStruct_in.Cgrid;
    % Adding name:
    names{lnames_in+1} = name_in;
    
    
    
    %% Checking grid is correct, modifying it otherwise:
    
    
    
    %% Filing new data in Tgrid and Cgrid matrices:
    
    
    
    %% Filing output structure:
    
    gridStruct_out = struct;
    gridStruct_out.method = method;
    gridStruct_out.names = names;
    gridStruct_out.xgrid = xgrid;
    gridStruct_out.ygrid = ygrid;
    gridStruct_out.zgrid = zgrid;
    gridStruct_out.Tgrid = Tgrid;
    gridStruct_out.Cgrid = Cgrid;
    
    
    


end
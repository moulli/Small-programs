function gridStruct_out = create_gridStruct(method, name_in, zcoord_in, cor_in, increment, totgrid)

%% Function that create gridStruct structure, based on a first dataset.
%
%  Data is provided through inputs cor_in and zcoord_in which represent
%  respectively correlation and coordinates. Then, grid is created to
%  include all the values and information is computed.
%
%
%% Inputs:
%
%  --method: name of the method employed to fill gridStruct.
%  --name_in: string, name of the new file added.
%  --zcoord_in: (n x 3) matrix, n being the number of neurons, and
%    coordinates should be taken from the ZBrain Atlas.
%  --cor_in: array of length n, with correlation coefficient for each
%    neuron provided.
%  --increment: can either be a number or a vector of length 3. Represents
%    the increment in the grid in the x, y and z dimensions.
%  --totgrid (optional): user can provide the algorithm with a grid already
%    computed. totgrid must be a cell of length 3, containing xgrid, ygrid
%    and zgrid. If just a vector is provided, program will assume xgrid, 
%    ygrid and zgrid are the same. 
%
%
%% Output:
%
%  --gridStruct_out: structure containing fields 'method', 'names', 'xgrid', 
%    'ygrid', 'zgrid', 'Tgrid' and 'Cgrid'. The first field is the method 
%    employed for the analysis, second field is the name of the HDF5 
%    already computed, third, fourth and fifth fields are information on 
%    the grid, the sixth one if a cell containing for each HDF5 all the 
%    neurons for each part of the grid, and seventh field is the mean
%    correlation for each part of the grid. Structure is created based on
%    given data.



    %% Initialization:
    
    % Indication:
    tic
    fprintf('\n\nLaunching function create_gridStruct.\n');
    % Checking coordinates and correlations:
    [nneu, n3d] = size(zcoord_in);
    if n3d ~= 3; error('Please provide coordinates with 3 columns for the 3 dimensions.'); end
    if ~isvector(cor_in); ('Please provide correlation as a vector.'); end
    if nneu ~= length(cor_in); error('Please provide coordinates and correlations with same size.'); end
    % Checking increment:
    if ~isvector(increment)
        error('Please provide increment as an array.')
    elseif sum(length(increment) == [1, 3]) == 0
        error('increment must be of length 1 or 3.')
    end
    increment = reshape(increment, 1, length(increment));
    increment = increment .* ones(1, 3);
    % Input grid:
    if nargin == 5
        tggiven = 0;
        totgrid = cell(1, 3);
    else
        if ~iscell(totgrid) || isvector(totgrid)
            totgrid = {totgrid};
        elseif ~iscell(totgrid) || ~isvector(totgrid)
            error('Please provide totgrid as a vector cell.')
        end
        if length(totgrid) == 1
            totgrid = {totgrid{1}, totgrid{1}, totgrid{1}};
        elseif length(totgrid) ~= 3
            error('Please provide totgrid with the right format.')
        end
        totgrid = reshape(totgrid, 1, 3);
        tggiven = 1;
    end  
    % Creating structure:
    gridStruct_out = struct;
    gridStruct_out.method = string(method);
    gridStruct_out.names = {name_in};
    
    
    
    %% Computing grid fields of output structure:
    
    % Indication:
    fprintf('Initialization completed in %.2f seconds.\n', toc);
    % Computing extrema:
    extrema = [min(zcoord_in); max(zcoord_in)];
    % Building grid:
    lgrid = zeros(1, 3);
    for i = 1:3 % for-loop on [x, y, z]
        if tggiven == 0
            totgrid{i} = extrema(1, i):increment(i):extrema(2, i);
        end
        lgrid(i) = length(totgrid{i});
    end
    % Assigning corresponding fields:
    xgrid = totgrid{1};
    gridStruct_out.xgrid = xgrid;
    ygrid = totgrid{2};
    gridStruct_out.ygrid = ygrid;
    zgrid = totgrid{3};
    gridStruct_out.zgrid = zgrid;
    
    
    
    %% Computing Tgrid and Cgrid based on input data:
    
    % Indication:
    fprintf('Grid computed in %.2f seconds, launching loop.\n', toc);
    % Building cell and matrix:
    Tgrid = cell(lgrid(1)-1, lgrid(2)-1, lgrid(3)-1, 1);
    Cgrid = zeros(lgrid(1)-1, lgrid(2)-1, lgrid(3)-1, 1);
    % Filling Tgrid and Cgrid based on grid:
    for ix = 1:(lgrid(1)-1)
        for iy = 1:(lgrid(2)-1)
            for iz = 1:(lgrid(3)-1)
                xtemp = (xgrid(ix) <= zcoord_in(:, 1) & zcoord_in(:, 1) < xgrid(ix+1));
                ytemp = (ygrid(iy) <= zcoord_in(:, 2) & zcoord_in(:, 2) < ygrid(iy+1));
                ztemp = (zgrid(iz) <= zcoord_in(:, 3) & zcoord_in(:, 3) < zgrid(iz+1));
                Ttemp = find(xtemp & ytemp & ztemp);
                Tgrid{ix, iy, iz, 1} = Ttemp;
                cortemp = mean(cor_in(Ttemp));
                if isnan(cortemp); cortemp = 0; end
                Cgrid(ix, iy, iz, 1) = mean(cortemp);   
            end
        end
        % Indication:
        if mod(ix, 4) == 0
            fprintf('For-loop %.2f %% completed in %.2f seconds.\n', [100*ix/(lgrid(1)-1), toc]);
        end
    end
    % Assinging corresponding fields:
    gridStruct_out.Tgrid = Tgrid;
    gridStruct_out.Cgrid = Cgrid;
    % Indication:
    fprintf('Function create_gridStruct ended in %.2f seconds.\n', toc);
    
    
end
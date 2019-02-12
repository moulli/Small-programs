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
    
    % Indication:
    tic
    fprintf('\n\nLaunching function add_to_gridStruct.\n');
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
    lTgrid_in = size(gridStruct_in.Tgrid, 4);
    lCgrid_in = size(gridStruct_in.Cgrid, 4);
    lx = length(gridStruct_in.xgrid);
    ly = length(gridStruct_in.ygrid);
    lz = length(gridStruct_in.zgrid);
    if lnames_in ~= lTgrid_in; error('Number of names does not correspond to number of examples in Tgrid field.'); end
    if lnames_in ~= lCgrid_in; error('Number of names does not correspond to number of examples in Cgrid field.'); end
    if lTgrid_in > 1
        if length(size(gridStruct_in.Tgrid)) ~= 4
            error('Tgrid field must be a 4 dimension matrix.')
        end
        if length(size(gridStruct_in.Cgrid)) ~= 4
            error('Cgrid field must be a 4 dimension matrix.')
        end
        if ~isequal([lx-1, ly-1, lz-1, lTgrid_in], size(gridStruct_in.Tgrid))
            error('Dimensions between grid and Tgrid cell do not match.')
        end
        if ~isequal([lx-1, ly-1, lz-1, lCgrid_in], size(gridStruct_in.Cgrid))
            error('Dimensions between grid and Cgrid matrix do not match.')
        end
    end
    % Checking coordinates and correlations:
    [nneu, n3d] = size(zcoord_in);
    if n3d ~= 3; error('Please provide coordinates with 3 columns for the 3 dimensions.'); end
    if ~isvector(cor_in); ('Please provide correlation as a vector.'); end
    if nneu ~= length(cor_in); error('Please provide coordinates and correlations with same size.'); end
    
    % Duplicating data from structure:
    method = gridStruct_in.method;
    names = gridStruct_in.names;
    xgrid = gridStruct_in.xgrid; xgrid = reshape(xgrid, 1, length(xgrid));
    ygrid = gridStruct_in.ygrid; ygrid = reshape(ygrid, 1, length(ygrid));
    zgrid = gridStruct_in.zgrid; zgrid = reshape(zgrid, 1, length(zgrid));
    Tgrid = gridStruct_in.Tgrid;
    Cgrid = gridStruct_in.Cgrid;
    % Getting increments:
    xincr = mean(gradient(xgrid));
    yincr = mean(gradient(ygrid));
    zincr = mean(gradient(zgrid));
    % Adding name:
    names{lnames_in+1} = name_in;
    
    
    
    %% Checking grid is correct, modifying it otherwise:
    
    % Indication:
    fprintf('Initialization completed in %.2f seconds.\n', toc);
    % Getting extremum values:
    extrema = [min(zcoord_in); max(zcoord_in)];
    % Checking grids and adding values if necessary:
    incr = [xincr, yincr, zincr];
    totgrid = {xgrid, ygrid, zgrid};
    for i = 1:3 % for-loop on [x, y, z]
        % Minimum values:
        if extrema(1, i) < totgrid{i}(1)
            mintemp = fliplr(totgrid{i}(1):-incr(i):extrema(1, i));
            mintemp = [mintemp(1)-incr(i), mintemp];
            lmt = length(mintemp);
            add_dim = zeros(1, 4);
            add_dim(i) = lmt;
            if lTgrid_in > 1
                Tgridtemp = cell(size(Tgrid) + add_dim);
                Cgridtemp = zeros(size(Tgrid) + add_dim);
            else
                Tgridtemp = cell([size(Tgrid), 1] + add_dim);
                Cgridtemp = zeros([size(Tgrid), 1] + add_dim);
            end
            Tgridtemp(add_dim(1)+1:end, add_dim(2)+1:end, add_dim(3)+1:end, :) = Tgrid;
            Tgrid = Tgridtemp;
            Cgridtemp(add_dim(1)+1:end, add_dim(2)+1:end, add_dim(3)+1:end, :) = Cgrid;
            Cgrid = Cgridtemp;
            totgrid{i} = [mintemp, totgrid{i}];
        end
        % Maximum values:
        if extrema(2, i) > totgrid{i}(end)
            maxtemp = totgrid{i}(end):incr(i):extrema(2, i);
            maxtemp = [maxtemp, maxtemp(end)+incr(i)];
            lmt = length(maxtemp);
            add_dim = zeros(1, 4);
            add_dim(i) = lmt;
            if lTgrid_in > 1
                Tgridtemp = cell(size(Tgrid) + add_dim);
                Cgridtemp = zeros(size(Tgrid) + add_dim);
            else
                Tgridtemp = cell([size(Tgrid), 1] + add_dim);
                Cgridtemp = zeros([size(Tgrid), 1] + add_dim);
            end
            Tgridtemp(1:size(Tgrid, 1), 1:size(Tgrid, 2), 1:size(Tgrid, 3), :) = Tgrid;
            Tgrid = Tgridtemp;
            Cgridtemp(1:size(Cgrid, 1), 1:size(Cgrid, 2), 1:size(Cgrid, 3), :) = Cgrid;
            Cgrid = Cgridtemp;
            totgrid{i} = [totgrid{i}, maxtemp];
        end
    end
    % Getting back xgrid, ygrid and zgrid:
    xgrid = totgrid{1};
    ygrid = totgrid{2};
    zgrid = totgrid{3};
        
        
%         % For xgrid:
%         if extrema(1, 1) < xgrid(1)
%             xminadd = fliplr(xgrid(1):-xincr:extrema(1, 1));
%             xminadd = [xminadd(1)-xincr, xminadd];
%             Tgrid = [cell(length(xminadd), ly, lz, lTgrid_in); Tgrid];
%             Cgrid = [zeros(length(xminadd), ly, lz, lCgrid_in); Cgrid];
%             xgrid = [xminadd, xgrid];
%             lx = length(xgrid);
%         end
%         if extrema(2, 1) > xgrid(end)
%             xmaxadd = fliplr(xgrid(end):xincr:extrema(2, 1));
%             xmaxadd = [xmaxadd, xmaxadd(end)+xincr];
%             Tgrid = [Tgrid; cell(length(xmaxadd), ly, lz, lTgrid_in)];
%             Cgrid = [Cgrid; zeros(length(xmaxadd), ly, lz, lCgrid_in)];
%             xgrid = [xgrid; xmaxadd];
%             lx = length(xgrid);
%         end
%         % For ygrid:
%         if extrema(1, 2) < ygrid(1)
%             yminadd = fliplr(ygrid(1):-yincr:extrema(1, 2));
%             yminadd = [yminadd(1)-yincr, yminadd];
%             Tgrid = [cell(lx, length(yminadd), lz, lTgrid_in); Tgrid];
%             Cgrid = [zeros(lx, length(yminadd), lz, lCgrid_in); Cgrid];
%             ygrid = [yminadd, ygrid];
%             ly = length(ygrid);
%         end
%         if extrema(2, 2) > ygrid(end)
%             ymaxadd = fliplr(ygrid(end):yincr:extrema(2, 2));
%             ymaxadd = [ymaxadd, ymaxadd(end)+yincr];
%             Tgrid = [Tgrid; cell(lx, length(yminadd), lz, lTgrid_in)];
%             Cgrid = [Cgrid; zeros(lx, length(yminadd), lz, lTgrid_in)];
%             ygrid = [ygrid; ymaxadd];
%             ly = length(ygrid);
%         end
%         % For zgrid:
%         if extrema(1, 3) < zgrid(1)
%             mintemp = fliplr(zgrid(1):-zincr:extrema(1, 3));
%             mintemp = [mintemp(1)-zincr, mintemp];
%             [xtemp, ytemp, ztemp, dtemp] = size(Tgrid);
%             Tgridtemp = cell(xtemp, ytemp, ztemp+length(mintemp), dtemp);
%             Tgridtemp(:, :, end-ztemp+1:end, :) = Tgrid;
%             Tgrid = Tgridtemp;
%             Cgridtemp = zeros(xtemp, ytemp, ztemp+length(mintemp), dtemp);
%             Cgridtemp(:, :, end-ztemp+1:end, :) = Cgrid;
%             Cgrid = Cgridtemp;
%             ygrid = [yminadd, ygrid];
%             ly = length(ygrid);
%         end
%         if extrema(2, 2) > ygrid(end)
%             ymaxadd = fliplr(ygrid(end):yincr:extrema(2, 2));
%             ymaxadd = [ymaxadd, ymaxadd(end)+yincr];
%             Tgrid = [Tgrid; cell(lx, length(yminadd), lz, lTgrid_in)];
%             Cgrid = [Cgrid; zeros(lx, length(yminadd), lz, lTgrid_in)];
%             ygrid = [ygrid; ymaxadd];
%             ly = length(ygrid);
%         end
    
    
    
    %% Filling new data to Tgrid and Cgrid matrices:
    
    % Indication:
    fprintf('Grid analysis completed in %.2f seconds, launching loop.\n', toc);
    % Assigning new cell for Tgrid, and new matrix for Cgrid, with new room:
    if lTgrid_in > 1
        Tgridtemp = cell(size(Tgrid) + [0, 0, 0, 1]);
        Cgridtemp = zeros(size(Cgrid) + [0, 0, 0, 1]);
    else
        Tgridtemp = cell([size(Tgrid), 1] + [0, 0, 0, 1]);
        Cgridtemp = zeros([size(Cgrid), 1] + [0, 0, 0, 1]);
    end
    Tgridtemp(:, :, :, 1:end-1) = Tgrid;
    Tgrid = Tgridtemp;
    Cgridtemp(:, :, :, 1:end-1) = Cgrid;
    Cgrid = Cgridtemp;
    % Algorithm:
    [lx, ly, lz, ~] = size(Tgrid);
    for ix = 1:lx
        for iy = 1:ly
            for iz = 1:lz
                xtemp = (xgrid(ix) <= zcoord_in(:, 1) & zcoord_in(:, 1) < xgrid(ix+1));
                ytemp = (ygrid(iy) <= zcoord_in(:, 2) & zcoord_in(:, 2) < ygrid(iy+1));
                ztemp = (zgrid(iz) <= zcoord_in(:, 3) & zcoord_in(:, 3) < zgrid(iz+1));
                Ttemp = find(xtemp & ytemp & ztemp);
                Tgrid{ix, iy, iz, end} = Ttemp;
                cortemp = mean(cor_in(Ttemp));
                if isnan(cortemp); cortemp = 0; end
                Cgrid(ix, iy, iz, end) = mean(cortemp);
            end
        end
        % Indication:
        if mod(ix, 4) == 0
            fprintf('For-loop %.2f %% completed in %.2f seconds.\n', [100*ix/lx, toc]);
        end
    end    
    
    
    
    %% Filing output structure:
    
    % Filling structure:
    gridStruct_out = struct;
    gridStruct_out.method = method;
    gridStruct_out.names = names;
    gridStruct_out.xgrid = xgrid;
    gridStruct_out.ygrid = ygrid;
    gridStruct_out.zgrid = zgrid;
    gridStruct_out.Tgrid = Tgrid;
    gridStruct_out.Cgrid = Cgrid;
    % Indication:
    fprintf('Function add_to_gridStruct ended in %.2f seconds.\n', toc);
    
    
end
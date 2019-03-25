classdef ZBraingrid < handle
    
%% Class to easily manipulate grid data.
%  
%  This class possesses the same properties as gridStruct had fields. Its
%  methods are creation, adding another example, plotting all or plotting
%  some of the examples.



    %% Properties:
    
    properties (SetAccess = private)
        
        %% Global information:
        % Method employed, simple string:
        method
        % Names of the examples, paths and comments:
        names
        paths
        comments
        
        %% Grid:
        % increment, size and (x, y, z)-grids:
        increment
        gridsize
        xgrid
        ygrid
        zgrid
        
        %% Data based on grid:
        % Correlations as a vector:
        Zcorvect
        % Classified data into grid, 4D cell:
        Zneurons
        % Correlations into grid, 4D matrix:
        Zcorrelations
        % Indication on number of neurons per example in grid, 4D matrix:
        Zneuron_number
        
    end
    
    
    
    %% Methods:
    
    % Dynamic:
    methods
        
        %% Constructor:
        function obj = ZBraingrid(method_in, increment_in) 
        % Constructor of class, fills method and creates grids.
            % Checking input:
            if nargin ~= 2
                error('Please provide method name and increment.')
            end
            if ~isvector(increment_in)
                error('Please provide increment as a vector.')
            elseif length(increment_in) == 1
                increment_in = increment_in .* ones(1, 3);
            end
            if length(increment_in) ~= 3
                error('Please provide increment as a scalar or as a 1x3 vector.')
            end
            % Adding method name:
            obj.method = string(method_in);
            % Building grid:
            obj.increment = increment_in;
            obj.xgrid = 0:increment_in:(0.5 + increment_in(1));
            obj.ygrid = 0:increment_in:(1.05 + increment_in(2));
            obj.zgrid = 0:increment_in:(0.3 + increment_in(3));
            % Building rest of properties:
            objsize = [length(obj.xgrid), length(obj.ygrid), length(obj.zgrid)];
            obj.gridsize = [objsize, 0];
            obj.Zcorvect = {};
            obj.Zneurons = cell(objsize - 1);
            obj.Zcorrelations = zeros(objsize - 1);
            obj.Zneuron_number = zeros(objsize - 1);   
            % Indication:
            % fprintf('ZBrain grid object created.\n');
        end
        
        %% Plus operation:
        function onew = plus(o1, o2)
        % Function that allows to add two ZBraingrids together.
            % Check methods and increments are the same:
            if ~isequal(string([o1.method]), string([o2.method]))
                error('The two objects should have the same methods.')
            elseif ~isequal([o1.increment], [o2.increment])
                error('The two objects should have the same increments.')
            end
            % Create new object:
            onew = ZBraingrid([o1.method], [o1.increment]);
            % Fill object with values from o1 and o2:
            onew.names = [[o1.names]; [o2.names]];
            onew.paths = [[o1.paths]; [o2.paths]];
            onew.comments = [[o1.comments]; [o2.comments]];
            onew.Zcorvect = [[o1.Zcorvect]; [o2.Zcorvect]];
            onew.Zneurons = cat(4, [o1.Zneurons], [o2.Zneurons]);
            onew.Zcorrelations = cat(4, [o1.Zcorrelations], [o2.Zcorrelations]);
            onew.Zneuron_number = cat(4, [o1.Zneuron_number], [o2.Zneuron_number]);
            onew.gridsize = size(onew.Zcorrelations);
        end
        
        %% Indexing operation:
        function onew = subsref(obj, Sin)
        % Function that allows indexing from a ZBraingrid object.
            switch Sin(1).type
                case '.'
                    if length(Sin) == 1
                        onew = builtin('subsref', obj, Sin);
                    elseif length(Sin) == 2
                        otemp = builtin('subsref', obj, Sin(1));
                        onew = builtin('subsref', otemp, Sin(2));
                    else
                        error('ZBraingrid:subsref',...
                              'Not a supported subscripted reference')
                    end
                case '()'
                    if length(Sin) == 1
                        % Build onew:
                        onew = ZBraingrid(obj.method, obj.increment);
                        onew.names = builtin('subsref', obj.names, Sin);
                        onew.paths = builtin('subsref', obj.paths, Sin);
                        onew.comments = builtin('subsref', obj.comments, Sin);
                        onew.Zcorvect = builtin('subsref', obj.Zcorvect, Sin);
                        Sin_b = struct('type', '()', 'subs', {[':', ':', ':', Sin.subs]});
                        onew.Zneurons = builtin('subsref', obj.Zneurons, Sin_b);
                        onew.Zcorrelations = builtin('subsref', obj.Zcorrelations, Sin_b);
                        onew.Zneuron_number = builtin('subsref', obj.Zneuron_number, Sin_b);
                        onew.gridsize = size(onew.Zcorrelations);
                        return
                    else
                        error('ZBraingrid:subsref',...
                              'Not a supported subscripted reference')
                    end
                case '{}'
                    error('ZBraingrid:subsref',...
                          'Not a supported subscripted reference')
            end
        end
        
        %% Absolut value:
        function onew = abs(obj, premean)
        % Function that takes the absolut value of the correlation.
        % Takes absolute value after mean across several examples.
        % If premean (optional) == 'premean', then absolute before mean:
        % this takes more into account neurons that react to absolute stim.
            % Checking premean:
            if nargin == 1
                premean = 0;
            elseif nargin == 2 && premean == "premean"
                premean = 1;
            elseif nargin ~= 2 || premean ~= "premean"
                error('Optional argument must be premean.')
            end
            % Building new object:
            onew = ZBraingrid(obj.method, obj.increment);
            onew.names = obj.names;
            onew.paths = obj.paths;
            onew.comments = obj.comments;
            onew.xgrid = obj.xgrid;
            onew.ygrid = obj.ygrid;
            onew.zgrid = obj.zgrid;
            onew.Zneurons = obj.Zneurons;
            onew.Zneuron_number = obj.Zneuron_number;
            % Filling positive correlations to leftover properties:
            [lx, ly, lz, ld] = size(obj.Zcorrelations);
            onew.Zcorvect = cell(ld, 1);
            onew.Zcorrelations = zeros(lx, ly, lz, ld);
            if premean == 0
                onew.Zcorvect = obj.Zcorvect;
                onew.Zcorrelations = builtin('abs', obj.Zcorrelations);
            else
                for id = 1:ld
                    ztemp = builtin('abs', obj.Zcorvect{id});
                    onew.Zcorvect{id} = ztemp;
                    for ix = 1:lx
                        for iy = 1:ly
                            for iz = 1:lz
                                cortemp = mean(ztemp(onew.Zneurons{ix, iy, iz, id}));
                                if isnan(cortemp); cortemp = 0; end
                                onew.Zcorrelations(ix, iy, iz, id) = cortemp;
                            end
                        end
                    end
                end
            end
            onew.gridsize = size(onew.Zcorrelations);
        end
        
        %% Length:
        function lnew = length(obj)
        % Function that provides length of ZBraingrid object.
            lnew = builtin('length', obj.Zcorvect);
        end
        
        %% Fill size for new object version:
        function fillsize(obj)
            obj.gridsize = size(obj.Zcorrelations);
        end
        
        %% Saving:
        function onew = saveobj(obj)
        % Function that compresses ZBraingrid object and saves it.
            % Creating new temporary object:
            onew = ZBraingrid(obj.method, obj.increment);
            onew.names = obj.names;
            onew.paths = obj.paths;
            onew.comments = obj.comments;
            onew.Zcorvect = obj.Zcorvect;
            % Find non-zero values:
            findind = find(obj.Zcorrelations ~= 0);
            lfind = length(findind);
            % Fill properties with lighter matrices:
            onew.gridsize = {obj.gridsize, findind}; % temporary 2 in 1
            onew.Zcorrelations = obj.Zcorrelations(findind);
            onew.Zneuron_number = obj.Zneuron_number(findind);
            % Replace cell with sparse matrix for Zneurons:
            laycell = obj.Zneurons(findind);
            maxcell = 0;
            for i = 1:lfind
                maxcell = max([maxcell, length(laycell{i})]);
            end
            neu_temp = zeros(lfind, maxcell);
            for i = 1:lfind
                neu_temp(i, 1:length(laycell{i})) = laycell{i}';
            end
            onew.Zneurons = sparse(neu_temp);
        end
        
        %% Add dataset to object:
        addDataset(obj, dataset_in);
        
        %% Plot all correlations, averaged over all datasets:
        plotAll(obj, varargin);
        
        %% Plot a subset of the correlations:
        plotSome(obj, subset, varargin);
        
        %% Flatten ZBraingrid object across all examples:
        onew = flatten(obj, opt_comment);
        
        %% Clean duplicates, if the same dataset is present more than once:
        cleanDuplicates(obj);
        
        %% Create a new object with lower increment:
        onew = downIncrement(obj, new_increment);
        
        %% Create subset object based on a comment keyword research:
        onew = choseSubset(obj, subset);
        
    end
    
    
    % Static:
    methods (Static)
        
        %% Getting rid of some low correlation neurons:
        [Vrid, varargout] = static_ridNeurons(Vin, ridvalues, varargin);
        
        %% Arranging color for plotting:
        Ccolor = static_corr2col(Ccorrelation, varargin)
        
        %% Loading:
        function onew = loadobj(loaded)
        % Function that decompresses ZBraingrid object from pathname.
            % Basic information:
            onew = ZBraingrid(loaded.method, loaded.increment);
            onew.names = loaded.names;
            onew.paths = loaded.paths;
            onew.comments = loaded.comments;
            onew.gridsize = loaded.gridsize{1};
            onew.Zcorvect = loaded.Zcorvect;
            % Computed information:
            sizetemp = loaded.gridsize{1};
            indtemp = loaded.gridsize{2};
            onew.Zcorrelations = zeros(sizetemp);
            onew.Zcorrelations(indtemp) = loaded.Zcorrelations;
            onew.Zneuron_number = zeros(sizetemp);
            onew.Zneuron_number(indtemp) = loaded.Zneuron_number;
            onew.Zneurons = cell(sizetemp);
                % Getting rid of the zeros:
                neutemp = num2cell(loaded.Zneurons);
                neutemp(find(loaded.Zneurons) == 0) = {[]}; % to have the same as initial empty values
                neutemp = num2cell(neutemp, 1);
                onew.Zneurons(indtemp) = strcat(neutemp{:});
        end
        
    end
    
    
end
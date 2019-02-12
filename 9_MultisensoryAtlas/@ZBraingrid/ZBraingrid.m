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
        % increment and (x, y, z)-grids:
        increment
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
            obj.Zcorvect = {};
            obj.Zneurons = cell(objsize - 1);
            obj.Zcorrelations = zeros(objsize - 1);
            obj.Zneuron_number = zeros(objsize - 1);   
            % Indication:
            fprintf('ZBrain grid object created.\n');
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
        end
        
        %% Indexing operation:
        function onew = subsref(obj, Sin)
        % Function that allows indexing from a ZBraingrid object.
            switch Sin(1).type
                case '.'
                    if length(Sin) == 1
                        onew = builtin('subsref', obj, Sin);
                    elseif length(Sin) == 2 || Sin(2).type == "()"
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
        
        %% Add dataset to object:
        addDataset(obj, dataset_in);
        
        %% Plot all correlations, averaged over all datasets:
        plotAll(obj, varargin);
        
        %% Plot a subset of the correlations:
        plotSome(obj, subset, varargin);
        
        %% Flatten ZBraingrid object across all examples:
        onew = flatten(obj, opt_comment);
        
        %% clean
        
        %% lower_increment
        
    end
    
    
    % Static:
    methods (Static)
        
        %% Getting rid of some low correlation neurons:
        [Vrid, varargout] = static_ridNeurons(Vin, ridvalues, varargin);
        
        %% Arranging color for plotting:
        Ccolor = static_corr2col(Ccorrelation, varargin)
        
    end
    
    
end
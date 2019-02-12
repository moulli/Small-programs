classdef ZBraingrid < handle
    
%% Class to easily manipulate grid data.
%  
%  This class possesses the same properties as gridStruct had fields. Its
%  methods are creation, adding another example, plotting all or plotting
%  some of the examples.



    %% Properties:
    
    properties
        
        % Method employed, simple string:
        method
        % Names of the examples, paths and comments:
        names
        paths
        comments
        
        % increment and (x, y, z)-grids:
        increment
        xgrid
        ygrid
        zgrid
        
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
        
        % Constructor:
        function obj = ZBraingrid(method_in, increment_in) 
        % Constructor of class, fills method and creates grids.
            % Checking input:
            if nargin ~= 2
                error('Please provide method name and increment.')
            end
            % Adding method name:
            obj.method = string(method_in);
            % Building grid:
            obj.increment = increment_in;
            obj.xgrid = 0:increment_in:(0.5 + increment_in);
            obj.ygrid = 0:increment_in:(1.05 + increment_in);
            obj.zgrid = 0:increment_in:(0.3 + increment_in);
            % Building rest of properties:
            objsize = [length(obj.xgrid), length(obj.ygrid), length(obj.zgrid)];
            obj.Zcorvect = {};
            obj.Zneurons = cell(objsize - 1);
            obj.Zcorrelations = zeros(objsize - 1);
            obj.Zneuron_number = zeros(objsize - 1);   
            % Indication:
            fprintf('ZBrain grid object created, with increment %f.\n', increment_in);
        end
        
        % Add dataset to object:
        addDataset(obj, dataset_in);
        
        % Plot all correlations, averaged over all datasets:
        plotAll(obj, varargin);
        
        % Plot some of the correlations:
        plotSome(obj);
        
    end
    
    
    % Static:
    methods (Static)
        
        % Getting rid of some low correlation neurons:
        [Vrid, varargout] = aux_ridNeurons(Vin, ridvalues, varargin);
        
        % Arranging color for plotting:
        Ccolor = aux_corr2col(Ccorrelation, varargin)
        
    end
    
    
end
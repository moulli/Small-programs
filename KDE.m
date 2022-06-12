classdef KDE < handle
    
%% Class to easily manipulate kernel density estimation
%  Created by Hippolyte on the 30/04/2020.
%  
%  The kernel density estimation is a probability function, continuous over
%  the space specified. It allows to have the probability that a specified
%  point belongs to the space defined by coordinates. In order to get this
%  probability distribution, we use a simple centered standard gaussian
%  kernel. The only parameter is h, the bandwidth, which is a 'smoothening'
%  parameters: the highest it is, the more space kde will cover. The lowest
%  it is the more discrete the final distribution will be.
%  See: https://en.wikipedia.org/wiki/Kernel_density_estimation
%
%  Properties:
%  - h [single/double]: bandwidth (see above)(mm).
%  - dim [integer]: dimension of data.
%  - len [integer]: number of points in the object.
%  - options [structure]: these are the plot options. They include the
%    isovalues ([0, 1]), the (grid) increment (mm), and the gridsize (mm).
%  - data [(n x dim) matrix]: points of interest (mm).
%  - kpdf [function]: continuous function returning probability density
%    value for each point in space associated to dim.
%
%  Methods:
%  - KDE: constructor. Takes bandwidth as parameter, whose default value is
%    0.0128mm (value found by Thijs).
%  - plus (+): when adding two KDE objects, this function checks that inner
%    bandwidth and data dimension is consistent, and if it is it creates a
%    new object with the data from the two added objects.
%  - bandwidth: changes the bandwidth.
%  - add: adds data to a KDE object.
%  - clear: deletes all data from KDE object.
%  - getpdf: defines PDF and stores it in the properties. In practice, this
%    function is not to be used by the user: new PDF is computed
%    automatically everytime there is a change in the parameters.
%  - pdf: given a set of points, computes the probability density for them
%    based on the PDF property.
%  - sample: without parameter, returns a point sampled from PDF, otherwise
%    returns as many sampled points as the input. In order to do that,
%    random numbers are picked, and indicate which data points we are going
%    to use. Then, random multivariate gaussian points are picked,
%    multiplied by the bandwidth and added to data points previously
%    picked.
%  - isovalues: changes the isovalues.
%  - increment: changes the increment.
%  - gridsize: changes the gridsize. IMPORTANT: for now, grid is defined
%    only on the positive part of space, i.e. it only goes from 0 to the
%    specified border in every direction. This could be changed if
%    necessary in the future.
%  - plot: plot using isovalues from options properties. I added an input
%    parser so that you can specify which isovalues, increment, and
%    gridsize you want without modifying it in the object.
%  - contour: contour using isovalues from options properties. I added an 
%    input parser so that you can specify which isovalues, increment, and
%    gridsize you want without modifying it in the object.
%  - project: get projection of the PDF from grid on a particular plan.


    %% Properties:
    
    % Private access
    properties (SetAccess = private)
        h % bandwidth
        dim % dimension of data
        len % number of points
        options % plot options
    end
    
    % Invisible
    properties (SetAccess = private, GetAccess = private)
        data % centers of data provided
        kpdf % continuous probability distribution function
    end 
    
    
    %% Methods
    
    % Dynamic
    methods
        
        %% Constructor
        function obj = KDE(h) 
        % Constructor of class. Defines bandwidth h.
            % Checking input:
            if nargin == 0
                h = 0.0128; % bandwidth value suggested by Thijs
            elseif ~isequal(size(h), [1, 1]) || h <= 0
                error('Bandwidth must be a positive number')
            end
            % add bandwidth property
            obj.h = h;
            obj.data = []; % for points of interest later
            obj.dim = 0;
            obj.len = 0;
            obj.getpdf();
            % plot options
            obj.options.isovalues = [0.4, 0.5, 0.7];
            obj.options.increment = 0.01 .* ones(1, 3);
            obj.options.gridsize = [0.496, 1.122, 0.276];
        end
        
        %% Add two KDE objects
        function onew = plus(o1, o2)
        % This function adds two KDE objects by adding their data properties
            if ~isequal(o1.h, o2.h)
                error('KDE objects do not have the same bandwidth')
            elseif ~isequal(o1.dim, o2.dim)
                error('Dimensions of datasets are not consistent')
            else
                onew = KDE(o1.h);
                onew.add(o1.data);
                onew.add(o2.data);
            end
        end        
        
        %% Change bandwidth
        function bandwidth(obj, h_new)
        % This function changes the bandwidth.
            if ~isequal(size(h_new), [1, 1]) || h_new <= 0
                error('Bandwidth must be a positive number')
            else
                obj.h = h_new;
                obj.getpdf();
            end
        end
        
        %% Add dataset to object
        function add(obj, coords)
        % This function adds coords to already existing data.   
            [l, d] = size(coords);
            if obj.dim == 0
                obj.data = coords;
                obj.dim = d;
                obj.len = obj.len + l;
                obj.getpdf();
            elseif d ~= obj.dim
                error('Dimension of new coordinates does not match')
            else
                obj.data = cat(1, obj.data, coords);
                obj.len = obj.len + l;
                obj.getpdf();
            end
        end
        
        %% Clear data
        function clear(obj)
        % This function clears obj.data    
            obj.data = [];
            obj.dim = 0;
            obj.len = 0;
            obj.getpdf();
        end
        
        %% Get PDF continuous function
        function getpdf(obj)
        % This function define obj.kpdf    
            % define kernel
            kernel = @(x) mean( mvnpdf( (x-obj.data) / obj.h ) / obj.h.^obj.dim );
            % function for all points
            obj.kpdf = @(x) cellfun( kernel, mat2cell( x, ones(size(x, 1), 1) ) );
        end
        
        %% Get probability distribution function
        function out = pdf(obj, pts)
        % This function computes the value of the PDF at points pts.   
            d = size(pts, 2);
            if obj.dim == 0
                error('No data to compute PDF from')
            elseif d ~= obj.dim
                error('Dimension of new coordinates does not match')
            else
                % define at given points
                out = obj.kpdf(pts);
            end
        end
        
        %% Sample from distribution
        function out = sample(obj, n)
        % This function returns n points sampled from pdf.
            if isempty(obj.data)
                error('No data to sample from')   
            else
                if nargin == 1
                    n = 1;
                end
                unisample = rand(n, 1);
                % get random gaussian from obj.data
                nc = size(obj.data, 1);
                gn = sum((0:nc-1)./nc <= unisample, 2);
                % get random normal value from these
                out = mvnrnd(zeros(1, obj.dim), eye(obj.dim), n) .* obj.h + obj.data(gn, :);  
            end
        end
        
        %% Modify isovalues
        function isovalues(obj, isovalues_new)
        % This function changes the isovalues in the plot options    
            [isox, isoy] = size(isovalues_new);
            if isox ~= 1 && isoy ~= 1
                error('Isovalues must be a number or an array')
            else
                obj.options.isovalues = reshape((isovalues_new), 1, length(isovalues_new));
            end
        end
        
        %% Modify grid increment
        function increment(obj, increment_new)
        % This function changes the isovalues in the plot options    
            [incx, incy] = size(increment_new);
            if incx ~= 1 && incy ~= 1
                error('Increment must be a number or an array')
            elseif incx == 1 && incy == 1
                obj.options.increment = increment_new .* ones(1, obj.dim);
            elseif incx == obj.dim || incy == obj.dim
                obj.options.increment = reshape((increment_new), 1, length(increment_new));
            else
                error('Dimension of increment must be 1 or dimension of data')
            end
        end
        
        %% Modify gridsize
        function gridsize(obj, gridsize_new)
        % This function changes the isovalues in the plot options    
            [gridx, gridy] = size(gridsize_new);
            if gridx ~= 1 && gridy ~= 1
                error('Gridsize must be a number or an array')
            elseif gridx == 1 && gridy == 1
                obj.options.gridsize = gridsize_new .* ones(1, obj.dim);
            elseif gridx == obj.dim || gridy == obj.dim
                obj.options.gridsize = reshape((gridsize_new), 1, length(gridsize_new));
            else
                error('Dimension of gridsize must be 1 or dimension of data')
            end
        end
        
        %% Plot isovalues
        function plot(obj, varargin)
        % This function plots isovalues based on options properties    
            % parse data
            p = inputParser;
            checksize = @(x) (size(x, 1) == 1 & size(x, 2) == 1) | (size(x, 1) == 1 & size(x, 2) == obj.dim) | (size(x, 1) == obj.dim & size(x, 2) == 1);
            isnotempty = @(x) ~isempty(x);
%             for i = 1:length(varargin)
%                 disp({varargin{i}, checksize(varargin{i}), isnotempty(varargin{i})})
%             end
            addOptional(p, 'increment', obj.options.increment, checksize);
            addOptional(p, 'isovalues', obj.options.isovalues, isnotempty);
            addOptional(p, 'gridsize', obj.options.gridsize, checksize);
            addOptional(p, 'Color', [1, 0, 0]);
            parse(p, varargin{:});
            % meshgrid
            x = 0:p.Results.increment:p.Results.gridsize(1);
            y = 0:p.Results.increment:p.Results.gridsize(2);
            z = 0:p.Results.increment:p.Results.gridsize(3);
            [X, Y, Z] = meshgrid(x, y, z);
            % define points
            gridpts = [X(:), Y(:), Z(:)];
            % evaluate
            prob = obj.pdf(gridpts);
            prob = reshape(prob, size(X));
            % divide by maximum value for easier manipulation with isovalues
            prob = prob ./ max(prob(:));
            % transparency for style
            alpha = linspace(0.1, 0.8, length(p.Results.isovalues));
            % for each isovalue
            for iv = 1:length(p.Results.isovalues)
                % get isovalue
                isov = p.Results.isovalues(iv);
                % get surface and patch
                surf = isosurface(X, Y, Z, prob, isov);
                pat = patch(surf);
                % define style
                isonormals(X, Y, Z, prob, pat);
                set(pat, 'FaceColor', p.Results.Color, 'EdgeColor', 'none', 'FaceAlpha', alpha(iv)); % set the color, mesh and transparency level of the surface
                if iv == 1
                    camlight(); 
                end
            end
        end
        
        %% Plot projection contour
        function contour(obj, plan, varargin)
        % This function plots isovalues based on options properties    
            % parse data
            p = inputParser;
            checksize = @(x) (size(x, 1) == 1 & size(x, 2) == 1) | (size(x, 1) == 1 & size(x, 2) == obj.dim) | (size(x, 1) == obj.dim & size(x, 2) == 1);
            isnotempty = @(x) ~isempty(x);
            addRequired(p, 'plan');
            addOptional(p, 'increment', obj.options.increment, checksize);
            addOptional(p, 'isovalues', obj.options.isovalues, isnotempty);
            addOptional(p, 'gridsize', obj.options.gridsize, checksize);
            parse(p, plan, varargin{:});
            % meshgrid
            x = 0:p.Results.increment:p.Results.gridsize(1);
            y = 0:p.Results.increment:p.Results.gridsize(2);
            z = 0:p.Results.increment:p.Results.gridsize(3);
            [X, Y, Z] = meshgrid(x, y, z);
            % define points
            gridpts = [X(:), Y(:), Z(:)];
            % evaluate
            prob = obj.pdf(gridpts);
            prob = reshape(prob, size(X));
            % divide by maximum value for easier manipulation with isovalues
            prob = prob ./ max(prob(:));
            % plot contour
            if p.Results.plan == 1
                prob = squeeze(max(prob, [], 2));
                xcontour = squeeze(Y(:, 1, :));
                ycontour = squeeze(Z(:, 1, :));
            elseif p.Results.plan == 2
                prob = squeeze(max(prob, [], 1));
                xcontour = squeeze(X(1, :, :));
                ycontour = squeeze(Z(1, :, :));
            elseif p.Results.plan == 3
                prob = squeeze(max(prob, [], p.Results.plan));
                xcontour = squeeze(X(:, :, 1));
                ycontour = squeeze(Y(:, :, 1));                
            end
            contour(xcontour, ycontour, prob, p.Results.isovalues)
        end
        
        %% Plot projection
        function out = project(obj, plan, varargin)
        % This function plots isovalues based on options properties    
            % parse data
            p = inputParser;
            checksize = @(x) (size(x, 1) == 1 & size(x, 2) == 1) | (size(x, 1) == 1 & size(x, 2) == obj.dim) | (size(x, 1) == obj.dim & size(x, 2) == 1);
            addRequired(p, 'plan');
            addOptional(p, 'increment', obj.options.increment, checksize);
            addOptional(p, 'gridsize', obj.options.gridsize, checksize);
            parse(p, plan, varargin{:});
            % meshgrid
            x = 0:p.Results.increment:p.Results.gridsize(1);
            y = 0:p.Results.increment:p.Results.gridsize(2);
            z = 0:p.Results.increment:p.Results.gridsize(3);
            [X, Y, Z] = meshgrid(x, y, z);
            % define points
            gridpts = [X(:), Y(:), Z(:)];
            % evaluate
            prob = obj.pdf(gridpts);
            prob = reshape(prob, size(X));
            % plot contour
            if p.Results.plan == 1
                out = squeeze(max(prob, [], 2));
            elseif p.Results.plan == 2
                out = squeeze(max(prob, [], 1));
            elseif p.Results.plan == 3
                out = squeeze(max(prob, [], p.Results.plan));               
            end
        end
        
    end
    
    
end
function plotKPDF(kpdf, varargin)

%% plotKPDF returns plot of isosurfaces that result from kpdf
%
%  First a grid is built, using increment value, and size provided. Then
%  kernel probabilities are computed at each grid point. Finally,
%  isosurfaces are displayed using matlab builtin function.
%
%  Inputs:
%  - kpdf [function]: associates to each point of space its probability to
%    belong to the kernel distribution.
%  - isovalues [array, optional]: all isovalues requested to plot
%    isosurfaces.
%  - increment [number, optional]: distance between two grid points (mm).
%  - gridsize [array, optional]: size of grid (mm).
%
%  Outputs:
%  - plot output.


    %% Check inputs
    
    % Default values
    defaultIsovalues = [0.4, 0.5, 0.7];
    defaultIncrement = 0.01;
    defaultGridsize = [0.496, 1.122, 0.276];
    
    % Input parser
    p = inputParser;
    addRequired(p, 'kpdf');
    addOptional(p, 'isovalues', defaultIsovalues);
    addOptional(p, 'increment', defaultIncrement);
    addOptional(p, 'gridsize', defaultGridsize);
    parse(p, kpdf, varargin{:});


    %% Build grid
    
    % meshgrid
    x = 0:p.Results.increment:p.Results.gridsize(1);
    y = 0:p.Results.increment:p.Results.gridsize(2);
    z = 0:p.Results.increment:p.Results.gridsize(3);
    [X, Y, Z] = meshgrid(x, y, z);
    
    
    %% Evaluate grid at each point
    
    % define points
    gridpts = [X(:), Y(:), Z(:)];
    
    % evaluate
    prob = kpdf(gridpts);
    prob = reshape(prob, size(X));
    
    % divide by maximum value for easier manipulation with isovalues
    prob = prob ./ max(prob(:));
    
    
    %% Plot isosurfaces
    
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
        set(pat, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', alpha(iv)); % set the color, mesh and transparency level of the surface
        if iv == 1
            camlight(); 
        end
    end
    
    
end
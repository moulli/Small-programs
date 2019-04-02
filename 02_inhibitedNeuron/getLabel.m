function out = getLabel(centroids, coordinates, labels, plotting, MaskDatabaseOutlines)
%% Function that will get the label associated to centroid coordinates.
%  The function will get the point that is the closest to each centroid in
%  coordinates. Then it will get the label associated through labels.
%  Optionally, the user can ask for a plot of all the regions associated to
%  the centroids. 'out' returns  a two column matrix, with two arrays 
%  comprising the label associated to each centroid, for major regions 
%  (1, 94, 114, 260 & 275) and minor regions.



    %% Get parameters:
    
    regions = 294;
    [icen, jcen] = size(centroids);
    [icoo, jcoo] = size(coordinates);
    [ilab, jlab] = size(labels);
    if icoo ~= ilab
        error('Please provide coordinates and labels with same number of rows')
    elseif jcen ~= 3
        error('Please provide centroids with 3 columns')
    elseif jcoo ~= 3
        error('Please provide coordinates with 3 columns')
    elseif jlab ~= regions
        error('Please provide labels with as many columns as there are regions')
    elseif ~(isequal(plotting, 'ON') || isequal(plotting, 'OFF'))
        error('Please provide plotting information as "ON" or "OFF"')  
    elseif isequal(plotting, 'ON') && nargin == 4
        error('Please provide MaskDatabaseOutlines.mat to have a plot of the regions')
    end
    
    
    
    %% Find closest point to each centroid:
    
    % Build objects to avoid for loops:
    dist3d = coordinates .* ones(icoo, jcoo, icen);
    permc = permute(centroids, [3, 2, 1]);
    % Compute distance and find minimum index:
    dtot = sum((dist3d - permc) .^ 2, 2);
    permd = permute(dtot, [1, 3, 2]);
    [~, indmin] = min(permd);
    
    
    
    %% Return the two labels for each centroid:
    
    out = cell(icen, 1);
    firstout = zeros(icen, 1);
    for i = 1:icen
        outemp = find(labels(indmin(i), :) == 1);
        if isempty(outemp) == 1
            continue
        else
            loutemp = length(outemp);
            sizeclus = zeros(loutemp, 1);
            for j = 1:loutemp
                masktemp = full(MaskDatabaseOutlines(:, outemp(j)));
                sizeclus(j) = sum(masktemp == 1);
            end
            [~, sortout] = sort(sizeclus);
            out{i} = outemp(sortout)';
            firstout(i) = outemp(sortout(1));
        end
    end
    
    
    
    %% Optional plot of outlines:
    
    % Meshgrid:
    width = 621; height = 1406; Zs = 138;
    x = 1:width;
    y = height:-1:1;
    z = 1:Zs;
    [X, Y, Z] = meshgrid(x, y, z);
    % Plotting:
    if plotting == 'ON'
        figure
        hold on
        grid on
        axis equal
        nnzfo = firstout(firstout~=0);
        unifo = length(unique(nnzfo));
        for i = 1:unifo
            founique = unique(nnzfo);
            masktemp = full(MaskDatabaseOutlines(:, founique(i)));
            sparsekeep = find(masktemp == 1);
            scatter3(X(sparsekeep), Y(sparsekeep), Z(sparsekeep), [], '.')
        end
    end  
    
    

end
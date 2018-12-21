function out = whichRegion(neurons, minmaxN, minmaxE, sparseData)
%% Function that will assign to each neuron in 'neurons' a brain region,
%  based on the data from MaskDatabase. Assignment is done based on the 
%  lowest distance between the neuron and a point in each region. The
%  lowest of these distances is the region associated to the neuron. If
%  several regions are picked, they will all be returned.
%  -'neurons' simply is the coordinates of the neurons we want to sort.
%  -'minmaxN' is a 2x3 matrix with the minimum and the maximum value of the
%   coordinates in all 3 directions.
%  -'minmaxE' is the same as 'minmaxN', but adapts to MaskDatabase.
%  -'sparseData' is actually MaskDatabase.mat from the zbrain atlas.
%  -'out' is given as a matrix, the length of 'neurons'. For each neuron, the
%   corresponding cell comprises all the minimum distances related to all 
%   the regions.

    

    %% Parameters:
    width = 621;
    height = 1406;
    Zs = 138;
    stot = width * height * Zs;
    nreg = 294;
    
    
    %% Informations from input data:
    [ineu, jneu] = size(neurons);
    [immN, jmmN] = size(minmaxN);
    [idata, jdata] = size(sparseData);
    if jneu ~= 3
        error('Please provide data as 3d points')
    elseif immN ~= 2 || jmmN ~= 3 || sum(minmaxN(1, :)<minmaxN(2, :)) ~= jmmN
        error('Please provide minimum & maximum values as a 2x3 matrix, with first row as minimum')
    elseif idata ~= stot || jdata ~= nreg
        error('Please provide MaskDatabase.mat as second argument')
    end
    
    
    %% Constructing out cell:
    out = zeros(ineu, nreg);
    
    
    %% Building appropriate and usable coordinates set for 'neurons':
    endsize = [width, height, Zs];
    neuCo = (neurons - minmaxN(1, :)) .* (minmaxE(2, :) - minmaxE(1, :)) ./ (minmaxN(2, :) - minmaxN(1, :)) + minmaxE(1, :);
    
    
    %% Constructing a meshgrid to convert the data to our convention:
    x = 1:width;
    y = height:-1:1;
    z = 1:Zs;
    [X, Y, Z] = meshgrid(x, y, z);
    
    
    %% Main algorithm - loop:
    for i = 1:ineu
        
        % Distance with each region's closest neuron:
        dist = zeros(1, nreg); 
        % Point considered:
        neuCoi = neuCo(i, :);
        
        % For each region:
        for j = 1:nreg
            % Taking associated region from MaskDatabase:
            regj = full(sparseData(:, j));
            % Finding location of region:
            locj = find(regj == 1);
            % Taking the subscript equivalent:
%             [X, Y, Z] = ind2sub(endsize, locj);
            % Reorganizing data to have the same convention:
            Ds = [X(locj), Y(locj), Z(locj)];
            % Computing the squared distance:
            distj = sum((Ds - neuCoi) .^ 2, 2);
            % Adding distance to distance arrow:
            dist(j) = min(distj);
            % Printing progress:
            fprintf('%i out of %i region, for %i out of %i neuron \n', [j, nreg, i, ineu]);
        end
        
        % Inserting distances into output:
        out(i, :) = dist;
        
    end

end
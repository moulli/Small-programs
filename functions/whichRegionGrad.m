function [out, distot] = whichRegionGrad(neurons, minmaxN, sparseData)
%% Function that will assign to each neuron in 'neurons' a brain region,
%  based on the data from MaskDatabase. Assignment is done based on the 
%  lowest distance between the neuron and a point in each region. The
%  lowest of these distances is the region associated to the neuron. If
%  several regions are picked, they will all be returned.
%  -'neurons' simply is the coordinates of the neurons we want to sort.
%  -'minmaxN' is a 2x3 matrix with the minimum and the maximum value of the
%   coordinates in all 3 directions.
%  -'sparseData' is actually MaskDatabase.mat from the zbrain atlas.
%  -'out' is given as a cell, the length of 'neurons'. For each neuron, the
%   corresponding cell comprises the assigned regions under the form of an
%   array.
%  We use a gradient technique to get to the lowest distance, to avoid
%  computing distances for all 120491388 points in MaskDatabase.

    

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
    out = cell(ineu, 1);
    
    
    %% Building appropriate and usable coordinates set for 'neurons':
    endsize = [width, height, Zs];
    neuCo = (neurons - minmaxN(1, :)) .* (endsize - 1) ./ (minmaxN(2, :) - minmaxN(1, :)) + 1;
    
    
    %% Main algorithm - loop:
    distot = cell(nreg, 1);
    for i = 1:ineu
        
        % Distance with each region's closest neuron:
        dist = zeros(nreg, 1); 
        
        % For each region:
        for j = 1:nreg
            % Taking associated region from MaskDatabase:
            regj = full(sparseData(:, j));
            % Finding location of region:
            locj = find(regj == 1);
            % Selecting random points:
            XR1 = locj(randperm(length(locj), 1));
            [XR1x, XR1y, XR1z] = ind2sub(endsize, XR1);
            distX = sqrt([XR1x, XR1y, XR1z] * neuCo(i, :)');
            distot{j} = [distot{j}; distX];
            XR2 = -1;
            % Gradient search:
            while nnz(XR1 ~= XR2) ~= 0 % if the 2 points are the same: stop
                % Look for the neighbours:
                [XR1x, XR1y, XR1z] = ind2sub(endsize, XR1);
                neigh = [XR1x + [-1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1; -1; 0; 1], ...
                         XR1y + [-1; -1; -1; 0; 0; 0; 1; 1; 1; -1; -1; -1; 0; 0; 0; 1; 1; 1; -1; -1; -1; 0; 0; 0; 1; 1; 1], ...
                         XR1z + [-1; -1; -1; -1; -1; -1; -1; -1; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1]];
                neigh = neigh(all(neigh, 2) & neigh(:, 1) <= width & neigh(:, 2) <= height & neigh(:, 3) <= Zs, :);
                neighind = sub2ind(endsize, neigh(:, 1), neigh(:, 2), neigh(:, 3));
                neigh = neigh(regj(neighind) == 1, :);
                if isempty(neigh) == 1
                    dist(j) = distX;
                    break
                else
                    distemp = sqrt(neigh * neuCo(i, :)');
                    [distXtemp, indist] = min(distemp);
                    distot{j} = [distot{j}; distXtemp];
                    if distXtemp == distX
                        break
                    else
                        % Updating:
                        distX = distXtemp;
                        XR2 = XR1;
                        XR1 = sub2ind(endsize, neigh(indist, 1), neigh(indist, 2), neigh(indist, 3));
                    end
                end
            end
            % Putting distance in distance array:
            dist(j) = distX;
            % Printing progress:
            fprintf('%i out of %i region, for %i out of %i neuron \n', [j, nreg, i, ineu]);
        end
        
        % Inserting closest region into output:
        out{i} = find(dist == min(dist));
        
    end

end
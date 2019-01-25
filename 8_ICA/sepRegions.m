function bzones = sepRegions(h5labels, h5coords, plotting)

%% Function that separates the different parts in the brain.
%
%  There are 10 different regions that are discriminated through this
%  algorithm. For each of the 5 following regions, there is left and right.
%  Telencephalon, diencephalon, mesencephalon, rhobencephalon and spinal
%  cord. Coordinates have to be given in a RAS orientation.
%
%
%% Parameters:
%
%  --h5labels: matrix with 294 columns, corresponding to 294 regions given
%    by the zbrain atlas.
%  --h5coords: matrix with 3 columns, corresponding to 3 dimensions for the
%    neurons of the brain.
%  --plotting: 'Plot' or 'plot' to plot zones, nothing otherwise.
%
%
%% Output:
%
%  --bzones: (5 x 2) cell with regions in rows, and left/right in columns.
%    For each cell of bzones, there is an array of indexes containing
%    neurons for each zone.



    %% Initialization:
    
    % Input:
    if nargin == 2; plotting = 'nope'; end
    [ml, nl] = size(h5labels);
    [mc, nc] = size(h5coords);
    if ml ~= mc; error('Please provide same number of rows for labels and coordinates.'); end
    if nl ~= 294; error('Please provide labels with 294 columns, like zbrain atlas.'); end
    if nc ~= 3; error('Please provide coordinates with 3 columns, for 3-dimension.'); end
    % Vector of zones that interest us:
    reglab = [275, 1, 94, 114, 260];
    % Output:
    bzones = cell(5, 2);
    
    
    %% Getting major label for each neuron, and solving problems:
    
    % Taking labels for each neuron:
    labelsnum = permute(h5labels .* (1:size(h5labels, 2)), [2, 3, 1]);
    labelslab = permute(sum(labelsnum == reglab), [3, 2, 1]);
    % Finding neurons misassigned:
    problems = find(sum(labelslab, 2) ~= 1);
    noproblems = find(sum(labelslab, 2) == 1);
    lprob = length(problems);
    labelslab(problems, :) = zeros(lprob, size(labelslab, 2));
    % K-nearest neighbours to find new labels:
    knn = 20;
    coordstemp = h5coords(noproblems, :);
    dist = sum((coordstemp - permute(h5coords(problems, :), [3, 2, 1])).^2, 2);
    [~, mindist] = sort(squeeze(dist));
    for i = 1:lprob
        labelstemp = sum(labelslab(noproblems(mindist(1:knn+1, i)), :));
        [~, maxlab] = max(labelstemp);
        labelslab(problems(i), maxlab) = 1;
    end
    
    
    %% Gathering neurons in a cell:
    
    % Computing discriminatory value for left/right separation:
    mid = mean(h5coords(:, 1));
    % Assigning values in bzones:
    for i = 1:5
        for j = 1:2
            switch j
                case 1
                    ktemp = ((labelslab(:, i) == 1) & (h5coords(:, 1) <= mid));
                case 2
                    ktemp = ((labelslab(:, i) == 1) & (h5coords(:, 1) > mid));
            end
            bzones{i, j} = find(ktemp == 1);
        end
    end
    
    
    %% Plotting if required:
    
    if lower(plotting) == 'plot'
        figure
        hold on
        for i = 1:5
            for j = 1:2
                scatter3(h5coords(bzones{i, j}, 1), h5coords(bzones{i, j}, 2), h5coords(bzones{i, j}, 3), [], rand(1, 3), '.')
                axis equal
            end
        end
    end



end
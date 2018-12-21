function out = dbscan(data, Params)

%% Density-based spatial clustering of applications w/ noise for 3D points.
%
%  This function clusters the data file, which must be provided as a set of
%  n-dimensions points. The clustering is done by keeping points whose 
%  entourage has at least a certain number of points. This gets rid of the 
%  points that are considered noise. Basic algorithm is improved with a new 
%  parameter: user can chose to add a weight to each point, to give more 
%  importance to certain points, and modulate this weight. Another
%  improvement is the fact that instead of having a rigid maximum distance,
%  the scalar given by user will define the standard deviation in the
%  smoothing gaussian.
%
%
%% Parameters:
%
%  --data: matrix of nD coordinates. Matrix must be provided with n columns,
%    and the number of rows is the number of points to analyse.
%  --Params: structure containing parameters for the function. Fields are:
%       --Params.dmax: maximum distance to compute the entourage of a neuron.
%       --Params.minpt: minimum number of points needed in the entourage of
%         a neuron so that is belongs to a cluster.
%       --Params.weight (default: 1 per each point): weights attibuted to
%         each point. Must be an array the same size as the number of rows
%         of the data matrix.
%       --Params.pweight (default: 1): influence of Params.weight.
%       --Params.plot (default: 0): if 2D or 3D data, plotting if 
%         Params.plot == 1. 
%       --Params.fprint (default: 5000): print progress after Params.fprint 
%         iterations done. 
%       --Params.monitor (default: 0): gives information on the algorithm.
%
%
%% Output:
%
%  --out: vector, same length as number of points in data (ie same length
%    as number of rows of data). Contains -1 if corresponding point has 
%    been labeled as noise, and the number of the cluster it has been 
%    assigned to otherwise.



    %% Initialization
    
    % Indication:
    tic
    fprintf('\n\nStarting program dbscan, for clustering data. \n');  
    % Input parameters:
    [ndata, mdata] = size(data);
    if ~isfield(Params, 'weight')
        Params.weight = ones(ndata, 1);
    else 
        [wx, wy] = size(Params.weight);
        if wy > wx; Params.weight = Params.weight'; [wx, wy] = size(Params.weight); end
        if wx ~= ndata || wy ~= 1
            error('Please provide Params.weight as an array with as length the number of points')
        end
    end
    if ~isfield(Params, 'pweight'); Params.pweight = 1; end
    if ~isfield(Params, 'plot') || (Params.plot == 1 && sum(mdata == [2, 3]) == 0)
        Params.plot = 0; end
    if ~isfield(Params, 'fprint'); Params.fprint = 5000; end
    if ~isfield(Params, 'monitor'); Params.monitor = 0; end
    if ~isscalar(Params.dmax)
        error('Please provide maximum distance as a scalar')
    elseif ~isscalar(Params.minpt)
        error('Please provide minimum number of points as a scalar')
    end
    % Normalizing weights:
    if length(unique(Params.weight)) ~= 1
        Params.weight = 4 * Params.pweight * (Params.weight - min(Params.weight)) ./ (max(Params.weight) - min(Params.weight)) - 2;
    end
    % Initialize variables and output:
    data = [(1:ndata)', data, Params.weight]; % to keep track
    datai = data;
    curclust = zeros(ndata, 1);
    nextinline = []; % which points are still to be analysed
    out = zeros(ndata, 1);
    
    
    
    %% DBscan clustering:
    
    % Indication:
    fprintf('\nInitialization done in %.3f seconds, starting loop for each point. \n\n', toc);
    iteration = 0;
    % Main loop:
    while nnz(out) < ndata
        % Monitoring how the algorithm works:
        if Params.monitor == 1
            unout = length(unique(out(out ~= 0 & out ~= -1)));
            if isempty(unout) == 1; unout = 0; end
            monitoring = [iteration, unout, nnz(out), nnz(curclust), length(nextinline), size(datai, 1)]
        end
        % Chosing point:
        if ~isempty(nextinline)
            % New point in the cluster:
            pti = nextinline(1);
            nextinline = nextinline(2:end);
        else
            if nnz(curclust) ~= 0
                % Updating clusters:
                out(curclust == 1) = max(out) + 1; 
                if isempty(find(out == 0, 1))
                    break
                end
                % Updating data:
                comparison = sum(datai(:, 1) == find(curclust == 0)', 2);
                datai = datai(comparison == 1, :);
                curclust = zeros(ndata, 1);
            end
            % Picking new point:
            restpti = find(out == 0);
            pti = data(restpti(randperm(length(restpti), 1)), 1);
        end
        % Computing distance:
        disttot = [datai(:, 1), sqrt(sum((datai(:, 2:end-1) - data(pti, 2:end-1)).^2, 2)), zeros(size(datai(:, 1)))]; % with index for facility
        disttot(:, 3) = datai(:, end) .* normpdf(disttot(:, 2), 0, Params.dmax) ./ normpdf(Params.dmax, 0, Params.dmax);
        dist = disttot(disttot(:, 2) < Params.dmax, :);
        lendist = sum(disttot(:, 3)) - normpdf(0, 0, Params.dmax) ./ normpdf(Params.dmax, 0, Params.dmax);
        % Assigning clusters:
        if lendist < Params.minpt && curclust(pti) == 0
            % Noise point:
            out(pti) = -1;
            datai = datai(datai(:, 1) ~= pti, :);
        elseif lendist >= Params.minpt
            if curclust(pti) == 0
                % Creating new cluster:
                curclust(pti) = 1;
            end
            % Adding points to the cluster:
            toadd = dist(:, 1);
            if isempty(nextinline) 
                nextinline = toadd(~sum(toadd == find(curclust ~= 0)', 2));
            else
                nextinline = [nextinline; toadd(~sum(toadd == nextinline', 2) & ~sum(toadd == find(curclust ~= 0)', 2))];
            end
            curclust(toadd) = 1;
        end
        % Indication:
        iteration = iteration + 1;
        if mod(iteration, Params.fprint) == 0
            fprintf('Iteration %.0f for %.0f points, done in %.3f seconds. \n', [iteration, ndata, toc]);
        end
    end
    
    
    
    %% Plotting if possible:
    
    if Params.plot == 1
        figure
        hold on
        axis equal
        grid on
        title('DBSCAN clustering on set of data', 'Interpreter', 'latex')
        xlabel('x-axis', 'Interpreter', 'latex')
        ylabel('y-axis', 'Interpreter', 'latex')
        clustmax = max(out);
        switch mdata
            case 2
                for i = 1:clustmax
                    scatter(data(out == i, 2), data(out == i, 3), ...
                             50, 'filled', 'MarkerFaceAlpha', 5/8, 'MarkerFaceColor', rand(1, 3))
                end
            case 3
                view(0, 90)
                zlabel('z-axis', 'Interpreter', 'latex')
                for i = 1:clustmax
                    scatter3(data(out == i, 2), data(out == i, 3), data(out == i, 4), ...
                             50, 'filled', 'MarkerFaceAlpha', 5/8, 'MarkerFaceColor', rand(1, 3))
                end
        end
    end
    
    % End of program indication:
    fprintf('\nFunction dbscan ended in %.3f seconds. \n', toc);
    


end
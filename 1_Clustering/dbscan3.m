function clusters = dbscan3(data, dmax, minpt, covariance, wcov, ~)
%% Density-based spatial clustering of applications with noise for 3D points
% - data must be provided as a matrix with 3 columns
% - dmax is the maximum distance from a point for the analysis
% - minpt is the number of points in the neighbourhood below which point is
% noise
% - covariance between signal and stimulus
% - wcov the weight we give to covariance (0: normal analysis)
% - plotting ('ON' or nothing) is whether or not we want to see algorithm
% live
% - clusters gives a column vector with the cluster associated to each point

    %% Clustering algorithm:
    
    % Plotting or not:
    if nargin == 5
        infoplot = 0;
    elseif nargin == 6
        infoplot = 1;
    end
    
    % Initializing:
    [ndata, mdata] = size(data);
    if mdata ~= 3
        error('Please provide data as a nx3 matrix')
    end
    data = [(1:ndata)', data, covariance]; % to keep track
    datai = data;
    clusters = zeros(ndata, 1);
    curclust = zeros(ndata, 1);
    nextinline = []; % which points are still to be analysed
    
    % Colors we will use:
    colmin = 0.1;
    colmax = 0.9;
    coloruse = [linspace(colmax, colmax, ceil(ndata/6))', linspace(colmin, colmax, ceil(ndata/6))', linspace(colmin, colmin, ceil(ndata/6))';
                linspace(colmax, colmin, ceil(ndata/6))', linspace(colmax, colmax, ceil(ndata/6))', linspace(colmin, colmin, ceil(ndata/6))';
                linspace(colmin, colmin, ceil(ndata/6))', linspace(colmax, colmax, ceil(ndata/6))', linspace(colmin, colmax, ceil(ndata/6))';
                linspace(colmin, colmin, ceil(ndata/6))', linspace(colmax, colmin, ceil(ndata/6))', linspace(colmax, colmax, ceil(ndata/6))';
                linspace(colmin, colmax, ceil(ndata/6))', linspace(colmin, colmin, ceil(ndata/6))', linspace(colmax, colmax, ceil(ndata/6))';
                linspace(colmax, colmax, ceil(ndata/6))', linspace(colmin, colmin, ceil(ndata/6))', linspace(colmax, colmin, ceil(ndata/6))'];
    
           
    % Algorithm:
    while nnz(clusters) < ndata
        % Chosing point:
        if ~ isempty(nextinline)
            % New point in the cluster:
            pti = nextinline(1);
            nextinline = nextinline(2:end);
        else
            if nnz(curclust) ~= 0
                % Updating clusters:
                clusters(curclust == 1) = max(clusters) + 1; 
                if isempty(find(clusters == 0, 1))
                    break
                end
                % Updating data:
                comparison = sum(datai(:, 1) == find(curclust == 0)', 2);
                datai = datai(comparison == 1, :);
                curclust = zeros(ndata, 1);
            end
            % Picking new point:
            restpti = find(clusters == 0);
            pti = data(restpti(randperm(length(restpti), 1)), 1);
        end
        % Computing distance:
        disttot = [datai(:, 1), sqrt(sum((datai(:, 2:4) - data(pti, 2:4)).^2, 2))]; % with index for facility
        disttot(find(disttot(:, 1) == pti), 2) = max(disttot(:, 2) + 2*dmax);
        dist = disttot(disttot(:, 2) < dmax, :);
        distcov = 1 + wcov .* (abs(covariance(dist(:, 1))) - median(abs(covariance))) ./ (max(abs(covariance)) - min(abs(covariance)));
        lendist = sum(distcov);
        % Assigning clusters:
        if lendist < minpt && curclust(pti) == 0
            % Noise point:
            clusters(pti) = -1;
            datai = datai(datai(:, 1) ~= pti, :);
        elseif lendist >= minpt
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
%         % Monitoring how the algorithm works:
%         [nnz(clusters), nnz(curclust), length(nextinline), sum(dist(:, 1) == 0), size(datai, 1)]
    end
    
    if infoplot == 1
        figure
        hold on
        axis equal
        view(0, 90)
        grid on
        title('DBSCAN clustering on 3D data', 'Interpreter', 'latex')
        xlabel('x-axis', 'Interpreter', 'latex')
        ylabel('y-axis', 'Interpreter', 'latex')
        zlabel('z-axis', 'Interpreter', 'latex')
        clustmax = max(clusters);
        for i = 1:clustmax
            scatter3(data(clusters == i, 2), data(clusters == i, 3), data(clusters == i, 4), ...
                     50, 'filled', 'MarkerFaceAlpha', 5/8, 'MarkerFaceColor', coloruse(randperm(size(coloruse, 1), 1), :))
        end
    end

end
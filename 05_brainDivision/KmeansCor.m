function [cent, maxclust, cent0] = KmeansCor(data, ncent, nite, display, init)
%% Functions that performs K-means algorithm on a set of data.
%  IT MAXIMIZES CORRELATION INSTEAD OF MINIMIZING EUCLIDIAN DISTANCE.
%  - data must be provided with examples as rows and features as columns.
%  - ncent: number of centroids we want.
%  - nite: number of iterations before returning results.
%  - display: 1 if we want progress printed, 0 otherwise.
%  - init[optional]: initial coordinates of centroids.
%  - cent: final coordinates of centroids, as a matrix.
%  - minclust: centroid associated to each point, as an arrow.
%  - cent0: initialization of centroids.



   %% Get data from parameters:
   
   [id, jd] = size(data);
   
   
   
   %% K-means initialization:
   
   % First clusters:
   tic % to have duration of algorithm
   switch nargin
       case 4
           cent = data(randperm(id, ncent), :);
       case 5
           [iin, jin] = size(init);
           if iin ~= ncent || jin ~= jd
               error('Please provide initialization with the right size')
           end
           cent = init;
   end
   cent0 = cent;
   % Assignment of points to cluster with maximum correlation:
   try
       clust3d = permute(cent, [3, 2, 1]);
       dist = permute(sum(data.*clust3d, 2) ./ sqrt(sum(data.^2, 2).*sum(clust3d.^2, 2)), [1, 3, 2])
   catch
       fprintf('\nHeavy files; kmeans will be in for loop (longer) \n');
       dist = zeros(id, ncent);
       for c = 1:ncent
           dist(:, c) = sum(data.*cent(c, :), 2) ./ sqrt(sum(data.^2, 2).*sum(cent(c, :).^2, 2));
       end
   end        
   [~, maxclust] = max(dist, [], 2);
   % Display time if required:
   if display == 1
       fprintf('Initialization done in %.3f seconds \n', toc);
   end
   
   
   
   %% Loop with new centroids and assignments:
   
   for i = 1:nite
       % New centroids:
       for c = 1:ncent
           cent(c, :) = mean(data(maxclust == c, :));
       end
       % New assignment of points to centroids:
       try
           clust3d = permute(cent, [3, 2, 1]);
           dist = permute(sum(data.*clust3d, 2) ./ sqrt(sum(data.^2, 2).*sum(clust3d.^2, 2)), [1, 3, 2]);
       catch
           dist = zeros(id, ncent);
           for c = 1:ncent
               dist(:, c) = sum(data.*cent(c, :), 2) ./ sqrt(sum(data.^2, 2).*sum(cent(c, :).^2, 2));
           end
       end
       [~, maxclust] = max(dist, [], 2);
       % Display number of iterations if required:
       if display == 1
           fprintf('Iteration %.0f out of %.0f in %.3f seconds \n', [i, nite, toc]);
       end
   end     
   


end
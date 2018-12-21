function [cent, minclust, cent0, rmsei] = iKmeans(data, ncent, nite, display, rmsestop)
%% Functions that performs K-means algorithm on a set of data.
%  Initialization is made through k++, and algo is completed by Expectation
%  Maximization algorithm.
%  - data must be provided with examples as rows and features as columns.
%  - ncent: number of centroids we want.
%  - nite: number of iterations before returning results.
%  - display: 1 if we want progress printed, 0 otherwise.
%  - rmsestop[optional]: stops algorithm if RMSE difference is lower.
%  - cent: final coordinates of centroids, as a matrix.
%  - minclust: centroid associated to each point, as an arrow.
%  - cent0: initialization of centroids.
%  - rmsei: RMSE for initialization and iterations.



   %% Get data from parameters:
   
   [id, jd] = size(data);
   tic % to have duration of algorithm
   % Display start of algorithm if required:
   if display == 1
       fprintf('\n\n\n\nK-means algorithm with k++ initialization launched \n\n');
   end
   
   
   
   %% K-means initialization:
   
   % First clusters with k++ algorithm:
   cent = zeros(ncent, jd);
   % k++ algorithm:
   cent(1, :) = data(randperm(id, 1), :);
   for i = 2:ncent
       distmin = zeros(id, i-1);
       for j = 1:(i-1)
           distmin(:, j) = sum((data - cent(j, :)).^2, 2);
       end
       distmin = min(distmin, [], 2).^2;
       distprob = [0; cumsum(distmin) ./ sum(distmin)];
       randnum = distprob - rand;
       randfind = find(randnum > 0);
       cent(i, :) = data(randfind(1)-1, :);
       % Display progress of initialization if required:
       if display == 1 && mod(i, 10) == 0
           fprintf('%.0f centroids initialized out of %.0f in %.3f seconds \n', [i, ncent, toc]);
       end
   end       
   cent0 = cent;
   % Assignment of points to closest cluster:
   try
       clust3d = permute(cent, [3, 2, 1]);
       dist = permute(sum((data - clust3d).^2, 2), [1, 3, 2]);
   catch
       fprintf('\nHeavy files; kmeans will be in for loop (longer) \n');
       dist = zeros(id, ncent);
       for c = 1:ncent
           dist(:, c) = sum((data - cent(c, :)).^2, 2);
       end
   end        
   [~, minclust] = min(dist, [], 2);
   % Compute RMSE:
   rmsei = zeros(nite+1, 1);
   rmsevect = zeros(id, jd);
   for c = 1:ncent
       rmsevect(minclust == c) = cent(c);
   end
   rmsei(1) = sum(sqrt(sum((data - rmsevect).^2, 2))) ./ id;
   % Display time if required:
   if display == 1
       fprintf('\nInitialization done in %.3f seconds, RMSE is %.3f \n\n', [toc, rmsei(1)]);
   end
   
   
   
   %% K-means loop with new centroids and assignments:
   
   for i = 1:nite
       % New centroids:
       for c = 1:ncent
           cent(c, :) = mean(data(minclust == c, :));
       end
       % New assignment of points to centroids:
       try
           clust3d = permute(cent, [3, 2, 1]);
           dist = permute(sum((data - clust3d).^2, 2), [1, 3, 2]);
       catch
           dist = zeros(id, ncent);
           for c = 1:ncent
               dist(:, c) = sum((data - cent(c, :)).^2, 2);
           end
       end
       [~, minclust] = min(dist, [], 2);
       rmsevect = zeros(id, jd);
       % Compute RMSE and stop if required:
       for c = 1:ncent
           rmsevect(minclust == c) = cent(c);
       end
       rmsei(i+1) = sum(sqrt(sum((data - rmsevect).^2, 2))) ./ id;
       if nargin == 5 && display == 1 && i >= 2 && abs(rmsei(i)-rmsei(i+1)) < rmsestop && abs(rmsei(i-1)-rmsei(i+1)) < rmsestop
           fprintf('\nRMSE is %.3f, algorithm stopped due to RMSE threshold \n\n', rmsei(i+1));
           rmsei = rmsei(all(rmsei, 2));
           break
       end
       % Display number of iterations if required:
       if display == 1
           fprintf('Iteration %.0f out of %.0f in %.3f seconds, RMSE is %.3f \n', [i, nite, toc, rmsei(i+1)]);
       end
   end    
   


end
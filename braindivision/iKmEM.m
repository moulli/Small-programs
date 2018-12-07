function [gparams, gsource, cent0, rmsei] = iKmEM(data, ncent, nite, display, rmsestop)
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
   
   
   
   %% EM loop based on found centroids:
   
   % Anounce beginning of EM if required:
   if display == 1
       fprintf('\n\nStarting EM algorithm \n\n');
   end
   % Initialization of gaussians:
   uniqueclust = unique(minclust);
   ncentf = length(uniqueclust);
   gparams = cell(ncentf, 3);
   % Creating vector in case determinant is numerically impossible:
   eigdet = zeros(ncentf, 1);
   for cn = 1:ncentf
       clustemp = (minclust == uniqueclust(cn));
       datatemp = data(clustemp, :);
       minclust(clustemp) = cn;
       gparams{cn, 1} = size(datatemp, 1) / id;
       gparams{cn, 2} = mean(datatemp);
       covtemp = cov(datatemp);
       gparams{cn, 3} = covtemp + 0.000001 * eye(jd);
       eigdet(cn) = mean(eig(gparams{cn, 3}));
   end
   % To have an estimate of distribution weight:
   eigdet = (min(eigdet) ./ eigdet) .^ jd;
   % Estimation of posterior:
   gsource = zeros(id, ncentf);
   for cn = 1:ncentf
       switch det(gparams{cn, 3})
           case 0
               gsource(:, cn) = eigdet(cn) * diag(exp(-0.5 * (data - gparams{cn, 2}) / gparams{cn, 3} * (data - gparams{cn, 2})'));
           otherwise
               gsource(:, cn) = gparams{cn, 1} .* mvnpdf(data, gparams{cn, 2}, gparams{cn, 3});
       end           
   end
   zerosumg = sum(gsource, 2) == 0; % In case some rows are empty.
   posterior = gsource ./ [sum(gsource, 2) + zerosumg];
   figure; image(posterior, 'CDataMapping', 'scaled')
   % Display time if required:
   if display == 1
       fprintf('\nInitialization of EM done (%.3f seconds) \n\n', toc);
   end
   % Starting loop:
   for i = 1:nite
       % Maximization step:
       for cn = 1:ncentf
           gparams{cn, 1} = mean(posterior(:, cn));
           postweight = posterior(:, cn) ./ sum(posterior(:, cn));
           gparams{cn, 2} = sum(postweight .* data);
           covtemp = (data - gparams{cn, 2})' * diag(postweight) * (data - gparams{cn, 2});
           gparams{cn, 3} = covtemp + 0.000001 * eye(jd);
           eigdet(cn) = mean(eig(gparams{cn, 3}));
       end
       % To have an estimate of distribution weight:
       eigdet = (min(eigdet) ./ eigdet) .^ jd;
       % Estimation step:
       for cn = 1:ncentf
           switch det(gparams{cn, 3})
               case 0
                   gsource(:, cn) = eigdet(cn) * diag(exp(-0.5 * (data - gparams{cn, 2}) / gparams{cn, 3} * (data - gparams{cn, 2})'));
               otherwise
                   gsource(:, cn) = gparams{cn, 1} .* mvnpdf(data, gparams{cn, 2}, gparams{cn, 3});
           end 
       end
       zerosumg = sum(gsource, 2) == 0; % In case some rows are empty.
       posterior = gsource ./ [sum(gsource, 2) + zerosumg];
      % Display time if required:
       if display == 1
           fprintf('Iteration %.0f out of %.0f in %.3f seconds for EM \n', [i, nite, toc]);
       end
   end
   


end
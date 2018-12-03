function out = corClust(nvalues, corval)
%% Function that clusters, based on covariance.
%  Creates groups of neurons, with similar activity.
%  - nvalues: matrix with parameters as columns, observations as rows.
%  - corval: value above which we consider 2 neurons to be linked.
%  - out: matrix containing neuron and group associated, in two rows.
    


    %% Initialization:
    
    % Time initialization:
    tic
    % out initialization:
    out = zeros(2, 0);
    % Matrix of covariance:
    covtemp = cov(nvalues');
    covtemp = covtemp - tril(covtemp);
    [mcheck, maxtemp] = max(covtemp(:));
    lcov = length(find(covtemp(:) >= corval));
    % To print progress:
    i = 0;
    
    
    
    %% Main algorithm:
    
    while mcheck >= corval
        % Printing progress:
        i = i+1;
        if mod(i, 10) == 0
            fprintf('Iteration %.0f, max correlation; %.3f, in %.2f seconds \n', [i, mcheck, toc]);
        end
        % Updating out:
        [itemp, jtemp] = ind2sub(size(covtemp), maxtemp);
        if sum(out(1, :) == itemp) == 1 && sum(out(1, :) == jtemp) == 1
            clust1 = out(2, out(1, :) == itemp);
            clust2 = out(2, out(1, :) == jtemp);
            [~, clustminind] = min([clust1, clust2]);
            switch clustminind
                case 1
                    out(2, out(2, :) == clust2) = clust1;
                    out(2, out(2, :) > clust2) = out(2, out(2, :) > clust2) - 1;
                case 2
                    out(2, out(2, :) == clust1) = clust2;
                    out(2, out(2, :) > clust1) = out(2, out(2, :) > clust1) - 1;
            end
        elseif sum(out(1, :) == itemp) == 1
            out = [out, [jtemp; out(2, out(1, :) == itemp)]];
        elseif sum(out(1, :) == jtemp) == 1
            out = [out, [itemp; out(2, out(1, :) == jtemp)]];
        else 
            out = [out, [itemp, jtemp; max([out(2, :), 0])+1, max([out(2, :), 0])+1]];
        end
        % Updating values:
        nvalues(itemp, :) = mean(nvalues([itemp, jtemp], :));
        nvalues(jtemp, :) = zeros(1, size(nvalues, 2));
        % Computing new covariance matrix:
        covtemp = cov(nvalues');
        covtemp = covtemp - tril(covtemp);
        [mcheck, maxtemp] = max(covtemp(:));
    end  
    
    

end
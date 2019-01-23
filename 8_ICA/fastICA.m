function [A, W, s] = fastICA(X, conv_or_nite, comp_or_var, funct)



    %% Dealing with input:
    
    % Indication:
    tic
    fprintf('\n\nFunction fastICA started\n\n');
    % Inputs:
    if nargin < 4; funct = 'tanh'; end
    if nargin < 3; comp_or_var = 1; end
    if nargin < 2; conv_or_nite = 100; end
    if floor(conv_or_nite) == conv_or_nite
        algo_ite = true;
        nite = conv_or_nite;
    else
        algo_ite = false;
        conv = conv_or_nite;
    end
    if floor(comp_or_var) == comp_or_var && comp_or_var ~= 1 && comp_or_var ~= 0
        algo_eig = true;
        component = comp_or_var;
    else
        algo_eig = false;
        varkeep = comp_or_var;
    end
    if comp_or_var < 0; error('Please provide retained variance or number of components comprised between 0 and 1, or positive integer'); end
    switch funct
        case 'tanh'
            g = @(u) tanh(u);
            gp = @(u) 1 - tanh(u).^2;
        case 'exp'
            g = @(u) u .* exp(-u.^2/2);
            gp = @(u) (1 - u.^2) .* exp(-u.^2/2);
        otherwise
            error('Please provide funct as tanh or exp.')
    end
    
    
    %% Centering data:    
    
    [~, n] = size(X);
    meanX = mean(X, 2);
    X = X - meanX;
    
    
    
    %% Whitening data:
    
    % Covariance matrix eigenvalues:
    covX = cov(X');
    [E, D] = eig(covX);
    % Determining which values to keep:
    [eigsort, eigind] = sort(diag(D), 'descend');
    if algo_eig
        eiglim = component;
    else
        eiglim = find((cumsum(eigsort) ./ sum(eigsort)) >= varkeep);
    end
    Emod = E(:, eigind(1:eiglim(1)));
    Dmod = diag(eigsort(1:eiglim(1)).^(-0.5));
    % New matrix:
    Xt = Dmod * Emod' * X;
    

    
    %% Main algorithm:
    
    % Indication:
    fprintf('\nCentering and whitening done in %.3f s; launching loop. \n', toc);
    % Initialization:
    A0 = randn(size(Xt, 1), size(Xt, 1));
    A0 = A0 / real((A0' * A0).^(0.5));
    % Algorithm:
    if algo_ite == true
        for i = 1:nite
            A0 = Xt * g(Xt' * A0) /n - mean(gp(Xt' * A0)) .* A0;
            A0 = A0 / real((A0' * A0)^(0.5));
            if mod(i, round(nite/10)) == 0
                fprintf('Iteration %.0f out of %.0f done in %.3f s. \n', [i, nite, toc]);
            end
        end
    else
        lastsave = 0;
        oldconv = 0;
        while true
            A1 = A0;
            A0 = Xt * g(Xt' * A0) /n - mean(gp(Xt' * A0)) .* A0;
            A0 = A0 / real((A0' * A0)^(0.5));
            newconv = mean(mean(abs(1 - (A0.*A1 ./ A0.^2))));
            if newconv <= conv || (newconv < oldconv && oldconv - newconv < 0.05)
                fprintf('Convergence criterion reached, or convergence evolution lowest than 0.05. \n');
                break
            elseif toc > lastsave+10
                fprintf('Convergence difference is %.6f for criterion %.6f done in %.3f s. \n', [newconv, conv, toc]);
                lastsave = lastsave+10;
            end
            oldconv = newconv;
        end
    end
    
    
    %% Dewhitening  and decentering:
    
    % First values:
    A = Emod / Dmod * A0;
    W = pinv(A);
    s = A \ (X + meanX);
    % 
    fprintf('\n\nFunction fastICA ended in %.0f s.\n\n', toc);


end
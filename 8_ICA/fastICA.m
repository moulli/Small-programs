function [A, W, s] = fastICA(X, conv_or_nite, varkeep, funct)



    %% Dealing with input:
    
    % Indication:
    tic
    fprintf('\n\nFunction fastICA started\n\n');
    % Inputs:
    if nargin < 4; funct = 'tanh'; end
    if nargin < 3; varkeep = 1; end
    if nargin < 2; conv_or_nite = 100; end
    if floor(conv_or_nite) == conv_or_nite
        algo_ite = true;
        nite = conv_or_nite;
    else
        algo_ite = false;
        conv = conv_or_nite;
    end
    if (varkeep > 1 || varkeep < 0); error('Please provide retained variance comprised between 0 and 1'); end
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
    
    [m, n] = size(X);
    meanX = mean(X, 2);
    X = X - meanX;
    
    
    
    %% Whitening data:
    
    % Covariance matrix eigenvalues:
    covX = cov(X');
    [E, D] = eig(covX);
    % Determining which values to keep:
    [eigsort, eigind] = sort(diag(D), 'descend');
    eiglim = find((cumsum(eigsort) ./ sum(eigsort)) >= varkeep);
    Emod = E(:, eigind(1:eiglim(1)));
    Dmod = diag(eigsort(1:eiglim(1)).^(-0.5));
    % New matrix:
    Xt = Dmod * Emod' * X;
    %Xt = Emod * Dmod * Emod' * X;
    

    
    %% Main algorithm:
    
    % Indication:
    fprintf('\nCentering and whitening done in %.3f s; launching loop. \n', toc);
    % Initialization:
    A0 = randn(m, size(Xt, 1));
    A0 = real((A0 * A0').^(-0.5)) * A0;
    %A0 = zeros(m);
    % Algorithm:
    if algo_ite == true
        for i = 1:nite
            A0 = Xt * g(Xt' * A0) /n - mean(gp(Xt' * A0)) * A0;
            %A0 = A0 ./ sqrt(sum(A0.^2));
            A0 = real((A0 * A0').^(-0.5)) * A0;
            if mod(i, round(nite/10)) == 0
                fprintf('Iteration %.0f out of %.0f done in %.3f s. \n', [i, nite, toc]);
            end
        end
        Aend = A0;
    else
        A1 = randn(m, size(Xt, 1)) * conv;
        A1 = real((A1 * A1').^(-0.5)) * A1;
        %A1 = randn(m, m) * conv;
        lastsave = 0;
        while sum(sum(abs(A1 - A0))) > conv
            A0 = A1;
            A1 = Xt * g(Xt' * A1) /n - mean(gp(Xt' * A1)) * A1;
            %A1 = A1 ./ sqrt(sum(A1.^2));
            A1 = real((A1 * A1').^(-0.5)) * A1;
            if toc > lastsave+10
                fprintf('Convergence difference is %.0f for criterion %.0f done in %.3f s. \n', [sum(sum(abs(A1 - A0))), conv, toc]);
                lastsave = lastsave+10;
            end
        end
        Aend = A1;
    end
    sum(isnan(Aend(:)))
    
    
    %% Dewhitening  and decentering:
    
    A = Emod / Dmod * Aend;
    %A = Emod * diag(diag(Dmod).^0.5) * Emod' * Aend;
    W = inv(A);
    s = A \ (X + meanX);
    fprintf('\n\nFunction fastICA ended in %.0f s.\n\n', toc);


end
function Fstat = computeFstat(signal, stimulus, time)
% Script explaining how to compute F-statistic for multilinear regression
% The F-statistic compares the multilinear regression to a degenerate
% regression, which means a regression using as the only regressor the
% constant signal equal to 1.


    %% Parameters
    
    % Degree of freedom linked to number of regressors (df1)
    p1 = 1; % dof for degenerate linear regression with only constant value as regressor
    p2 = 5; % dof for multilinear model used in neural analysis
    
    % Degree of freedom linked to number of points (df2)
    n = length(stimulus);
    
    
    %% Rearange signal and stimulus
    
    signal = reshape(signal, n, 1);
    stimulus = reshape(stimulus, n, 1);
    
    
    %% Building regressor for degenerate linear regression
    
    regressor1 = ones(n, 1);
    
    
    %% Building regressor for multilinear regression
    
    % Add regressors
    regressor2 = zeros(n, 5);
    regressor2(:, 1) = 1;
    regressor2(:, 2) = stimulus .* (stimulus > 0);
    regressor2(:, 3) = stimulus .* (stimulus < 0);
    regressor2(:, 4) = gradient(stimulus) .* (gradient(stimulus) > 0);
    regressor2(:, 5) = gradient(stimulus) .* (gradient(stimulus) < 0);
    regressor2(isnan(regressor2)) = 0;

    % Convolve with exponential kernel
    expdecay = 2.6;
    expkern =  exp(-(0:time:(20*expdecay))/expdecay);
    for i = 2:5
        regressorsi = conv(regressor2(:, i), expkern);
        regressor2(:, i) = regressorsi(1:size(regressor2, 1));
    end

    % Normalize
    for i = 2:5
%         % Classic normalization
%         regressors(:, i) = (regressors(:, i) - mean(regressors(:, i))) ./ std(regressors(:, i));
%         % Normalization using mode
%         regressors(:, i) = (regressors(:, i) - mode(regressors(:, i))) ./ std(regressors(:, i));
        % Normalization using mode and dividing by max value
        regressor2(:, i) = (regressor2(:, i) - mode(regressor2(:, i))) ./ max(abs(regressor2(:, i)));
    end
    
    
    %% Residual sum of squares for degenerate linear regression
    
    [coef1, ~, residual1] = regress(signal, regressor1);
    RSS1 = sum(residual1.^2);
    % NB: residual1 = signal - coef1.*regressor1
    
    
    %% Residual sum of squares for multilinear regression
    
    [coef2, ~, residual2] = regress(signal, regressor2);
    RSS2 = sum(residual2.^2);
    % NB: residual2 = signal - sum(coef2.*regressor2, 2)
    
    
    %% F-statistic value
    
    Fstat = ((RSS1-RSS2)/(p2-p1)) / (RSS2/(n-p2));
    


end
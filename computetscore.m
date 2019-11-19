function tscores = computetscore(signal, stimulus, time)
% Script explaining how to compute t-scores for multilinear 
% The t-scores are obtained for each coefficient by divinding it by its
% standard error. The standard error is itself computed as the square root
% of the residuals variance times the inverse of regressors'*regressors.
% The F-statistic compares the multilinear regression to a degenerate
% regression, which means a regression using as the only regressor the
% constant signal equal to 1.
    
    
    %% Rearange signal and stimulus
    
    n = length(stimulus);
    signal = reshape(signal, n, 1);
    stimulus = reshape(stimulus, n, 1);
    
    
    %% Building regressor for multilinear regression
    
    % Add regressors
    regressors = zeros(n, 5);
    regressors(:, 1) = 1;
    regressors(:, 2) = stimulus .* (stimulus > 0);
    regressors(:, 3) = stimulus .* (stimulus < 0);
    regressors(:, 4) = gradient(stimulus) .* (gradient(stimulus) > 0);
    regressors(:, 5) = gradient(stimulus) .* (gradient(stimulus) < 0);
    regressors(isnan(regressors)) = 0;

    % Convolve with exponential kernel
    expdecay = 2.6;
    expkern =  exp(-(0:time:(20*expdecay))/expdecay);
    for i = 2:5
        regressorsi = conv(regressors(:, i), expkern);
        regressors(:, i) = regressorsi(1:size(regressors, 1));
    end

    % Normalize
    for i = 2:5
        % Classic normalization
        regressors(:, i) = (regressors(:, i) - mean(regressors(:, i))) ./ std(regressors(:, i));
%         % Normalization using mode
%         regressors(:, i) = (regressors(:, i) - mode(regressors(:, i))) ./ std(regressors(:, i));
%         % Normalization using mode and dividing by max value
%         regressor2(:, i) = (regressor2(:, i) - mode(regressor2(:, i))) ./ max(abs(regressor2(:, i)));
    end
    
    
    %% Using normal equation to get coefficients
    
    X = regressors;
    y = signal;
    coefficients = (X' * X) \ X' * y;
    
    
    %% Computing residuals and residuals variance
    
    residuals = y - X * coefficients;
    var_residuals = var(residuals);
    
    
    %% Using variance to get standard error of coefficients
    
    coefficients_variance = var_residuals * pinv(X' * X);
    standard_error = sqrt(diag(coefficients_variance));
    
    
    %% Deducing t-score
    
    tscores = coefficients ./ standard_error;
    

end
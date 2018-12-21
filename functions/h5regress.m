function [out, variables, perinf] = h5regress(h5file, varin)

%% Function that performs linear regression on signals against variables.
%
%  With variables provided from HDF5 file, signal is averaged on the period 
%  of the variables if needed. Then a multilinear regression is performed 
%  against the provided variables.
%  NB: DATA ARE NORMALIZED BEFORE REGRESSION!
%
%
%% Parameters:
%
%  --h5file: HDF5 with a 'DFF' or 'Values' path, which is a matrix with as
%    many rows as there are observations, and as many columns as there are 
%    points per observation.
%  --varin: matrix against which to regress. Number of rows represent 
%    the number of variables, and the number of columns is the number of
%    points, which must be inferior to number of columns for signal.
%
%
%% Output:
%
%  --out: structure containing data on multilinear regression. Fields are:
%        --out.coef: coefficients for each neuron to variables.
%        --out.Fstat: 
%        --out.Pval:
%        --out.intercept:
%        --out.R2score:
%        --out.residual:
%        --out.samplevariance:
%  --variables: normalized input variables.
%  --perinf: matrix containing averaged normalized signal .


    
    %% Initialization:
    
    % Indication:
    tic
    fprintf('\n\nStarting program h5regress, for regression against variables. \n');
    % Getting files:
    try
        dff = h5read(h5file, '/Data/DFF');
    catch
        dff = h5read(h5file, '/Data/Values');
    end
    % Getting dimensions:
    [nneu, ntime] = size(dff);
    [nvarin, mvarin] = size(varin);
    % Normalizing data:
    varinn = (varin - mean(varin, 2)) ./ std(varin, [], 2);
    dffn = (dff - mean(dff, 2)) ./ std(dff, [], 2);
    
    
    %% Filling output:
    
    out = struct;
    out.coef = zeros(nneu, nvarin);
    out.Fstat = zeros(nneu, 1);
    out.Pval = zeros(nneu, 1);
    out.intercept = zeros(nneu, 1);
    out.R2score = zeros(nneu, 1);
    out.residual = zeros(nneu, mvarin);
    out.samplevariance = zeros(nneu, 1);
    
    
    %% Defining averaged signals:
    
    fullper = floor(ntime/mvarin);
    partper = mod(ntime, mvarin);
    dffntemp = zeros(nneu, mvarin);
    for j = 1:nneu
        meandffj = mean(reshape(dffn(j, 1:fullper*mvarin), mvarin, fullper), 2)';
        dffntemp(j, :) = fullper/(fullper+1)*meandffj + 1/(fullper+1)*[dffn(j, fullper*mvarin+1:end), meandffj(partper+1:end)];
    end
    dffn = dffntemp;
    % Mean signal:
    perinf = dffn;
    % Final variables:
    variables = [ones(size(varinn(1, :)')), varinn'];
    
    
    %% Multilinear regression:
    
    % Indication:
    fprintf('\nInitialization done in %.3f seconds, starting loop for each signal. \n\n', toc);
    % Main loop:
    for i = 1:nneu
        [b, ~, r, ~, stats] = regress(dffn(i, :)', variables);
        out.coef(i, :) = b(2:end)';
        out.Fstat(i) = stats(2);
        out.Pval(i) = stats(3);
        out.intercept(i) = b(1);
        out.R2score(i) = stats(1);
        out.residual(i, :) = r';
        out.samplevariance(i) = stats(4);
        % Indication:
        if mod(i, floor(nneu/10)) == 0
            fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nneu, toc]);
        end
    end      
    % End of program indication:
    fprintf('\nFunction h5regress ended in %.3f seconds. \n', toc);
    


end



function outkern = convInd(u, v, i1, i2)

%% Function that performs convolution and keep only indicated indices.
%
%
%% Parameters:
%  --u: first vector to convolve.
%  --v: second vector to convolve.
%  --i1: index to which we start to keep information.
%  --i2: index to which we stop to keep information.
%
%
%% Output:
%  --outkern: convolved vector, with i2-i1+1 components.



    %% Initialization:
    
    if ~isrow(u) && ~iscolumn(u) 
        error('Please provide first input as a vector')
    elseif ~isrow(v) && ~iscolumn(v)
        error('Please provide second input as a vector')
    elseif i1 < 1 || i1 > i2 || i2 > length(u)+length(v)-1
        error('Inputs must be as follow: 0 < input3 <= input4 < length(input1)+length(input2)')
    end
    
    
    %% Main code:
    
    outkern = conv(u, v);
    outkern = outkern(i1:i2);   



end
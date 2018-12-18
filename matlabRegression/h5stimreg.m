function [out, variables, perinf] = h5stimreg(h5file, Params)

%% Function that performs linear regression on signals against stimulus.
%
%  The stimulus is first averaged on a certain number of periods. Signals
%  are averaged the same way. Then a multilinear regression is performed
%  against positive part of stimulus, negative part of stimulus, positive
%  part of stimulus derivative, and negative part of stimulus derivative.
%  All of these variables are first convolved with an exponential kernel.
%  NB: DATA ARE NORMALIZED BEFORE REGRESSION!
%
%
%% Parameters:
%
%  --h5file: HDF5 file or link to HDF5 file. The file must have a
%    '/Data/Stimulus' path, a '/Data/Values' or '/Data/DFF' path, and a
%    '/Data/Times' path. Stimulus must be a (1 x time) array, DFF a 
%    (neuron x time) matrix, and Times a (1 x time) array.
%  --Params[optional]: structure containing parameters for the function. 
%    Fields are:
%       --Params.period (default: 0): if 0 taking all stimulus, else taking 
%         a certain number of periods.
%       --Params.nperiod (default: 1): number of periods on which to 
%         average.
%       --Params.expdecay (default: 2.6): time constant for exponential 
%         decay convolution.
%       --Params.fprint (default: 5000): print progress after Params.fprint 
%         iterations done.
%
%
%% Output:
%
%  --out: structure containing data on multilinear regression. Fields are:
%        --out.coef: coefficients for each neuron to 4 stimulus variables.
%        --out.Fstat: 
%        --out.Pval:
%        --out.intercept:
%        --out.R2score:
%        --out.residual:
%        --out.samplevariance:
%  --variables: positive stimulus, negative stimulus, positive stimulus
%    difference and negative stimulus difference convolved with the
%    exponential kernel, with a column of ones.
%  --perinf: cell with retained period, mean stimulus and mean signal.


    
    %% Initialization:
    
    % Indication:
    tic
    fprintf('\n\nStarting program h5stimreg, for regression against stimulus. \n');
    % Getting files:
    stim = h5read(h5file, '/Data/Stimulus');
    time = h5read(h5file, '/Data/Times');
    try
        dff = h5read(h5file, '/Data/DFF');
    catch
        dff = h5read(h5file, '/Data/Values');
    end
    % Getting dimensions:
    [nneu, ntime] = size(dff);
    if ntime ~= size(stim, 2) || ntime ~= size(time, 2)
        error('Please provide stimulus, times and signals with same number of time increments')
    end
    % Normalizing data:
    stimn = (stim - mean(stim)) ./ std(stim);
    dffn = (dff - mean(dff, 2)) ./ std(dff, [], 2);
    
    
    %% Filling structures:
    
    % If Params not provided, Params initialized here:
    if nargin == 1
        Params = struct;
    end
    % Input parameters:
    if ~isfield(Params, 'period'); Params.period = 0; end
    if ~isfield(Params, 'nperiod'); Params.nperiod = 1; end
    if ~isfield(Params, 'expdecay'); Params.expdecay = 2.6; end
    if ~isfield(Params, 'fprint'); Params.fprint = 5000; end
    % Defining output:
    out = struct;
    out.coef = zeros(nneu, 4);
    out.Fstat = zeros(nneu, 1);
    out.Pval = zeros(nneu, 1);
    out.intercept = zeros(nneu, 1);
    out.R2score = zeros(nneu, 1);
    out.residual = zeros(nneu, ntime);
    out.samplevariance = zeros(nneu, 1);
    
    
    %% Defining stimuli variables:
    
    % If Params.period == 1, computing mean period:
    perfin = ntime;
    if Params.period == 1
        % Finding period, using findpeaks:
        [~, indpeak] = findpeaks(stimn);
        indpeakb = indpeak(2:end) - indpeak(1);
        indpeakb(indpeakb > ntime/2) = [];
        lindpeakb = length(indpeakb);
        percor = zeros(lindpeakb, 1);
        pererr = zeros(lindpeakb, 1);
        for i = 1:lindpeakb
            stimn1 = stim(1:indpeakb(i));
            stimn2 = stim(indpeakb(i)+1:2*indpeakb(i));
            percor(i) = sum(stimn1.*stimn2) ./ sqrt(sum(stimn1.^2).*sum(stimn2.^2));
            pererr(i) = mean((stimn1-stimn2).^2); 
        end
        [~, indper] = max(percor ./ pererr);
        perfin = indpeakb(indper);
        % Taking number of periods given by Params.nperiod:
        if perfin * Params.nperiod <= ntime
            perfin = perfin .* Params.nperiod;
        end
        % Averaging stimulus and signal:
        fullper = floor(ntime/perfin);
        partper = mod(ntime, perfin);
        meanstim = mean(reshape(stimn(1:fullper*perfin), perfin, fullper), 2)';
        stimn = fullper/(fullper+1)*meanstim + 1/(fullper+1)*[stimn(fullper*perfin+1:end), meanstim(partper+1:end)];
        dffntemp = zeros(nneu, perfin);
        for j = 1:nneu
            meandffj = mean(reshape(dffn(j, 1:fullper*perfin), perfin, fullper), 2)';
            dffntemp(j, :) = fullper/(fullper+1)*meandffj + 1/(fullper+1)*[dffn(j, fullper*perfin+1:end), meandffj(partper+1:end)];
        end
        dffn = dffntemp;
        % New number of time increments:
        ntime = perfin;
        out.residual = zeros(nneu, ntime);
        % Indicating chosen period:
        fprintf('\nStimulus and signal averaged for period %.0f increments in %.3f seconds.', [perfin, toc]);
    end
    perinf = {perfin, stimn, dffn};
    % Taking positive and negative position and velocity:
    posStim = stimn .* (stimn > 0);
    negStim = stimn .* (stimn < 0);
    diff = gradient(stimn);
    diffn = (diff - mean(diff)) ./ std(diff);
    posDiff = diffn .* (diffn > 0);
    negDiff = diffn .* (diffn < 0);
    % Convolving with an exponential kernel:
    dt = mean(gradient(time));
    expkern =  exp(-(0:dt:(20*Params.expdecay))/Params.expdecay);
    pSKern = convInd(posStim, expkern, 1, ntime);
    nSKern = convInd(negStim, expkern, 1, ntime);
    pDKern = convInd(posDiff, expkern, 1, ntime);
    nDKern = convInd(negDiff, expkern, 1, ntime);
    % Final variables:
    variables = [ones(size(pSKern')), pSKern', nSKern', pDKern', nDKern'];
    
    
    %% Multilinear regression:
    
    % Indication:
    fprintf('\nInitialization done in %.3f seconds, starting loop for each neuron. \n\n', toc);
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
        if mod(i, 5000) == 0
            fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nneu, toc]);
        end
    end      
    % End of program indication:
    fprintf('\nFunction h5stimreg ended in %.3f seconds. \n', toc);
    


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
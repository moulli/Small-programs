function [data, info, N] = genSetConv(params)

%% Function that creates a set of data to train NN on deconvolution.
%
%  Based on the parameters given by the user, genSetConv provides a certain
%  number of data in order to be able to use them to train a neural
%  network. These data will be used to study deconvolution.
%
%
%% Parameters:
%
%  --params: structure containing parameters for the function. Fields are:
%       --params.nex: number of examples we want.
%       --params.pts: number of points per example we want.
%       --params.fs: sampling frequency.
%       --params.fr: maximum firing frequency.
%       --params.taur: rise time constant for convolution.
%       --params.taud: decay time constant for convolution.
%       --params.a: maximum height.
%       --params.b: baseline fluorescence for spikeless signal.
%       --params.noise: noise standart deviation.
%       --params.ssec (default: 1): number of spikes per second (> 0).
%       --params.per (default: 0): periodicity of spiking. Spiking
%         increases periodically.
%       --params.perinf (default: 0.1): influence of periodicity (must be
%         between 0 and 1).
%       --params.indic (default: 1000): a message is displayed for each of
%         the indic data created.
%  NB: {fs, fr, taur, taud, a, b, snr, ssec, per, perinf} CAN BE NUMBERS OR 
%  INTERVALS. IN THE CASE OF AN INTERVAL, RANDOM NUMBERS IN INTERVAL WILL 
%  BE TAKEN.
%
%
%% Output:
%
%  --data: (npts x nex) matrix containing the data.
%  --info: (1 x nex) matrix of structures with the info on fs, fr, taur, taud,
%    a, b, snr, ssec, per and perinf.


    
    %% Initialization:
    
    % Indication:
    tic
    fprintf('\n\nStarting program genSetConv, for creation of dataset. \n');
    % Getting input:
    if length(params.nex) ~= 1 || length(params.pts) ~= 1
        error('Please provide number of examples and number of points as a scalar')
    end
    if ~isfield(params, 'ssec'); params.ssec = 0.5; end
    if any(params.ssec < 0)
        error('Please provide ssec as a positive scalar')
    end
    if ~isfield(params, 'per'); params.per = 0; end
    if ~isfield(params, 'perinf'); params.perinf = 0.1; end
    if any(params.perinf) < 0 || any(params.perinf > 1)
        error('Please provide perinf as a scalar between 0 and 1')
    end
    % Defining labels to make it easier:
    labels = {'fs', 'fr', 'taur', 'taud', 'a', 'b', 'noise', 'ssec', 'per', 'perinf'};
    for i = 1:length(labels)
        ptemp = params.(labels{i});
        if length(ptemp) == 1
            params.(labels{i}) = ptemp * ones(1, 2);
        elseif length(ptemp) == 2
            params.(labels{i}) = sort(ptemp);
        else
            error('Please provide parameters as scalars or intervals between two points')
        end
    end
    if ~isfield(params, 'indic'); params.indic = 1000; end
    % Checking if taur < taud:
    tauc = params.taur - params.taud;
    if any(tauc(:) >= 0)
        error('Please provide a rise time constant inferior to the decay time constant')
    end
    % Building output:
    data = zeros(params.pts, params.nex);
    for i = 1:length(labels)
        infotemp = rand(1, params.nex);
        infotemp = infotemp * (params.(labels{i})(2) - params.(labels{i})(1)) + params.(labels{i})(1);
        info.(labels{i}) = infotemp;
    end
    % Indication:
    fprintf('\nInitialization done in %.3f seconds, starting data creation. \n\n', toc);
    
    
    
    %% Creation of data:

    % Exponential kernel function:
    expkern = @(t, tr, td) (exp(-t./td) - exp(-t./tr)) ./ max(exp(-t./td) - exp(-t./tr));
    % Now creating the dataset:
    for i = 1:params.nex
        % Getting parameters:
        fs = info.fs(i); fr = info.fr(i); taur = info.taur(i); taud = info.taud(i); 
        a = info.a(i); b = info.b(i); noise = info.noise(i); ssec = info.ssec(i); 
        per = info.per(i); perinf = info.perinf(i);
        % Constructing spiking train:
        xtemp = (0:(1/fr):((params.pts-1)/fs))';
        switch per
            case 0
                Ntemp = rand(size(xtemp));
            otherwise
                Ntemp = rand(size(xtemp)) + perinf * sin(2*pi * xtemp / per + rand * 2*pi);
        end
        Ntemp = double(Ntemp > exp(-ssec/fr));
        N(:, i) = Ntemp;
        % Constructing exponential kernel:
        expt = (0:(1/fr):(15*taud))';
        expk = expkern(expt, taur, taud);
        % Convolving:
        Ntemp = convInd(Ntemp, expk, 1, length(Ntemp));
        % Adding other parameters:
        Ntemp = a * Ntemp + noise * randn(size(Ntemp)) + b;
        % Finally sampling:
        xs = round(linspace(1, length(Ntemp), params.pts));
        data(:, i) = Ntemp(xs);
        % Indication:
        if mod(i, params.indic) == 0
            fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, params.nex, toc]);
        end
    end
    
    
    
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
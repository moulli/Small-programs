function createRegression(F)
% Computes regression with Focus input. If there are several stimuli, all
% stimuli available are regressed upon. Regression results are saved in the
% HDF5 folder.

    %% Recover information from Focus' HDF5
    
    % Recover HDF5 folder and path
    h5folder = F.dir('HDF5');
    h5path = fullfile(F.dir('HDF5'), [erase(F.name, ' ') '.h5']);
    
    % Recover neurons signals
    try 
        dff = h5read(h5path, '/Data/Brain/Analysis/DFFaligned');
    catch
        dff = h5read(h5path, '/Data/Brain/Analysis/DFF');
    end
    
    % Recover stimulus.i information
    h5infostim = h5info(h5path);
    for i = 1:size(h5infostim.Groups.Groups, 1)
        if h5infostim.Groups.Groups(i).Name == "/Data/Stimulus"
            break
        end
    end
    epath = cell(0, 1);
    for j1 = 1:size(h5infostim.Groups.Groups(i).Groups, 1)
        numpath = size(h5infostim.Groups.Groups(i).Groups(j1).Datasets, 1);
        for j2 = 1:numpath
            eptemp = fullfile(h5infostim.Groups.Groups(i).Groups(j1).Name, h5infostim.Groups.Groups(i).Groups(j1).Datasets(j2).Name);
            stimtemp = h5read(h5path, eptemp);
            if length(stimtemp) == size(dff, 2)
                epath = [epath; eptemp];
            end
        end
    end
    
    
    %% Compute regressions for each stimulus
    
    % Initialize output
    Regression = cell(0, length(epath));
    
    % Regression for each stimulus
    warning('off')
    for stimnum = 1:length(epath)
        % Load stimulus
        stim = h5read(h5path, epath{stimnum});
        stim = reshape(stim, length(stim), 1);
        % Compute regressors
        regressors = [ones(length(stim), 1), abs(expconv(stim.*(stim>0))), abs(expconv(stim.*(stim<0))), ...
                      abs(expconv(gradient(stim).*(gradient(stim)>0))), abs(expconv(gradient(stim).*(gradient(stim)<0)))];
        for i = 2:size(regressors, 2)
            regressors(:, i) = (regressors(:, i)-mean(regressors(:, i))) ./ std(regressors(:, i));
        end
        % Compute regression
        coefs = zeros(size(dff, 1), size(regressors, 2));
        R2 = zeros(size(dff, 1), 1);
        Fstat = zeros(size(dff, 1), 1);
        pvalue = zeros(size(dff, 1), 1);
        for i = 1:size(dff, 1)
            [coef, ~, ~, ~, stats] = regress(dff(i, :)', regressors);
            coefs(i, :) = coef';
            R2(i) = stats(1);
            Fstat(i) = stats(2);
            pvalue(i) = stats(3);
        end
        % Add output related to this stimulus
        Regression(stimnum).stimulus = epath{stimnum};
        Regression(stimnum).regressors = regressors;
        Regression(stimnum).coefs = coefs;
        Regression(stimnum).R2 = R2;
        Regression(stimnum).Fstat = Fstat;
        Regression(stimnum).pvalue = pvalue;
    end
    warning('on')
    
    
    %% Save output in HDF5 folder
    
    pathRegression = fullfile(h5folder, 'Regression.mat');
    save(pathRegression, 'Regression')
    
    
end
    


function out = expconv(signal, tau, dt)
% Function that convolves signal with exponential kernel.
    %% Default values:    
    if nargin == 2
        dt = 0.4;
    elseif nargin == 1
        tau = 2.6;
        dt = 0.4;
    elseif nargin ~= 3
        error('Please provide between 1 and 3 inputs.')
    end
    
    %% Building exponential kernel:
    texp = 0:dt:(20*tau);
    expkern = exp(-texp / tau);
       
    %% Convolving and keeping beginning of vector:    
    out = conv(signal, expkern);
    out = out(1:length(signal));    

end

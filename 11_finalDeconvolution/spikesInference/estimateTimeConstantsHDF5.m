function [tau_rises, tau_decays] = estimateTimeConstantsHDF5(h5_path, varargin)
% Use Tubiana's BSD algorithm to estimate GCaMP time constants.

    % Parallel BSD
    % ----------
    Palg = struct;
    Oalg = struct;
    Oalg.adaptive = 1;                  % Adaptative, refine given parameters
    Oalg.iterations = 5;                % Number of iterations for estimating parameters
    labelin = [];
    par = true;
    dffaligned = '/Data/Brain/Analysis/DFF';
    if nargin ~= 1
        lvar = length(varargin);
        if mod(lvar, 2) ~= 0
            error('Please provide inputs as pairs.')
        else
            inputs = ["tauRise", "tauDecay", "labels", "parfor", "adaptive", "iterations", "aligned"];
            for i = 1:2:lvar
                whichinput = find(varargin{i} == inputs);
                if isempty(whichinput)
                    error('Please use right keywords.')
                else
                    if whichinput == 1
                        Palg.tauRise = varargin{i+1};
                    elseif whichinput == 2
                        Palg.tauDecay = varargin{i+1};
                    elseif whichinput == 3
                        if ~isvector(varargin{i+1}) || ~isnumeric(varargin{i+1}) || max(varargin{i+1}) > 294 || min(varargin{i+1}) < 1 || mod(varargin{i+1}, 1) ~= 0
                            error('Please provide labels with right format.')
                        else
                            labelin = varargin{i+1};
                        end
                    elseif whichinput == 4
                        if sum(varargin{i+1} == [true, false]) == 0
                            error('Please provide parallel information as [true|false].')
                        else
                            par = varargin{i+1};
                        end
                    elseif whichinput == 5
                        if sum(varargin{i+1} == [false, true]) == 0
                            error('Adaptive parameter should be [false|true].')
                        else
                            Oalg.adaptive = varargin{i+1} * 1;
                        end
                    elseif whichinput == 6
                        if mod(varargin{i+1}, 1) ~= 0 || varargin{i+1} < 1
                            error('Number of iterations should be a positive integer.')
                        else
                            Oalg.iterations = varargin{i+1};
                        end
                    elseif whichinput == 7
                        if sum(varargin{i+1} == [false, true]) == 0
                            error('Choice for DF aligned should be [false|true].')
                        elseif varargin{i+1} == true
                            dffaligned = '/Data/Brain/Analysis/DFFaligned';
                        end
                    end
                end
            end
        end
    end
    fprintf('Initialization done.\n');
%     guess_rise = 0.5;
%     guess_decay = 3;       % educated guess to find time constants (seconds)
    
    % Main program
    % ----------
    
        % Getting DFF and time
        % ----------          
        dff = h5read(h5_path, dffaligned);  
        dff = dff';  
        dff = cast(dff, 'double'); 
        time = h5read(h5_path, '/Data/Brain/Times'); 
        time = cast(time, 'double'); 

        % Getting rid of too small and nan values
        % ----------
        dff(abs(dff) > 10*std(dff(:))) = 1e-3; 
        dff(isnan(dff)) = 1e-3; 
        dff(dff == 0) = 1e-3; 
    
        % Getting zones to keep
        % ----------
        if ~isempty(labelin)
            coord = h5read(h5_path, '/Data/Brain/ZBrainCoordinates');
            try
                labels = h5read(h5_path, '/Data/Brain/Labels');
            catch
                fprintf('Labels not found, computing labels.\n');
                labels = addLabels(coord);
                h5create(h5_path, '/Data/Brain/Labels', size(labels), 'Datatype', 'single');
                h5write(h5_path, '/Data/Brain/Labels', single(labels));
                h5writeatt(h5_path, '/Data/Brain/Labels', 'origin', 'ZBrain Atlas');
            end   
            % Recovering right label
            labels = labels(:, labelin);
            % Keeping just one side
            labels = (sum(labels, 2) > 0) & (coord(:, 1) < mean(coord(:, 1)));
            dff = dff(:, labels);
        end

        % Inference algorithm parameters
        % ----------                  % Struct of experimental conditions & decoding options.
        Oalg.Time = size(dff, 1);           % Number of time frames.
        Oalg.dt = mean(gradient(time));     % interval duration. (s)
        Oalg.nNeurons = size(dff, 2);       % Number of neurons.

        tic
        fprintf('Estimating time constants... \n');
        switch par
            case true
                [~, ~, ~, Pphys] = pBSD(dff, Oalg, Palg);
                delete(gcp('nocreate'));
            case false
                [~, ~, ~, Pphys] = BSD(dff, Oalg, Palg);
        end

        fprintf('Done (%2.2fs).\n', toc)

        tau_rises = Pphys.tauRise;
        tau_decays = Pphys.tauDecay;
        
    % Rest of code (optional)
    % ----------

%     figure;
%     subplot(121);
%     histogram(tau_rises, 100);
%     xlabel('\tau_{Rise} [s]');
%     subplot(122);
%     histogram(tau_decays, 100);
%     xlabel('\tau_{Decay} [s]');

    R.mean = mean(tau_rises);
    R.median = median(tau_rises);
    R.std = std(tau_rises);

    D.mean = mean(tau_decays);
    D.median = median(tau_decays);
    D.std = std(tau_decays);

    fprintf('Rising time : \n\t Mean = %12.8fs \n\t Median = %12.8fs \n\t Std = %12.8fs \n', ...
        R.mean, R.median, R.std);
    fprintf('Decay time : \n\t Mean = %12.8fs \n\t Median = %12.8fs \n\t Std = %12.8fs \n', ...
        D.mean, D.median, D.std);

end

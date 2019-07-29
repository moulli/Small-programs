function showProgress(iteration, maxvalue, division, barlength)
% Function that will fprintf a progression bar for any iterative algorithm.
%   -iteration: current iteration in the algorithm,
%   -maxvalue: maximum iteration in the algorithm,
%   -division (default 50): number of time progress bar will refresh,
%    division should be lower than maxvalue,
%   -barlength (default 50): length of progress bar.


    %% Possibility to change design:
    
    progress_sign = '>>>';
    progress_bar = '='; % this should be of length 1


    %% Default values:
    
    % Automatic values:
    if nargin == 3
        barlength = 50;
    elseif nargin == 2
        barlength = 50;
        division = min([50, maxvalue]);
    end
    % Check division:
    if division > maxvalue
        error('division should be lower than maxvalue.')
    % Check all values are integers:
    elseif floor(iteration)~= iteration
        error('Please provide iteration as an integer.')
    elseif floor(maxvalue) ~= maxvalue
        error('Please provide maxvalue as an integer.')
    elseif floor(division) ~= division
        error('Please provide division as an integer.')
    elseif floor(barlength)~= barlength
        error('Please provide barlength as an integer.')
    end

    
    %% Initial empty progress bar:
    
    % Resetting progress bar length:
    lsign = length(progress_sign);
    barlength_n = barlength - lsign + 1;
    % Initial progress bar:
    if iteration == 1
        fprintf('[');
        fprintf(progress_sign(end));
        fprintf(repmat(' ', 1, barlength_n-1));
        fprintf(']\n\n');
    end
    
    
    %% Iterative progress bar:
    
    % Define metric elements:
    roundtemp = round(maxvalue/division);
    floortemp = floor(iteration*barlength/maxvalue);
    % Conditional refresh:
    if iteration == maxvalue
        fprintf(repmat('\b', 1, barlength_n+5));
        fprintf('\n[');
        fprintf(repmat(progress_bar, 1, barlength_n));
        fprintf(']\n');
    elseif mod(iteration, roundtemp) == 0
        fprintf(repmat('\b', 1, barlength_n+3));
        fprintf(repmat(progress_bar, 1, floortemp-lsign+1));
        fprintf(progress_sign(max([1, lsign-floortemp]):end));
        fprintf(repmat('\b', 1, floortemp+1-barlength_n));
        fprintf(repmat(' ', 1, barlength_n-floortemp-1));
        fprintf(']\n\n');
    end


end

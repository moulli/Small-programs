function showProgress(iteration, maxvalue, division, barlength)
% Function that will fprintf a progression bar for any iterative algorithm.
%   -iteration: current iteration in the algorithm,
%   -maxvalue: maximum iteration in the algorithm,
%   -division (default 50): number of time progress bar will refresh,
%    division should be lower than maxvalue,
%   -barlength (default 50): length of progress bar.


    %% Default values:
    
    % Automatic values:
    if nargin == 3
        barlength = 50;
    elseif nargin == 2
        barlength = 50;
        division = 50;
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
    
    if iteration == 1
        fprintf('[>');
        fprintf(repmat(' ', 1, barlength-1));
        fprintf(']\n\n');
    end
    
    
    %% Iterative progress bar:
    
    % Define metric elements:
    roundtemp = round(maxvalue/division);
    floortemp = floor(iteration*barlength/maxvalue);
    % Conditional refresh:
    if iteration == maxvalue
        fprintf(repmat('\b', 1, barlength+5));
        fprintf('[');
        fprintf(repmat('=', 1, barlength));
        fprintf(']\n');
    elseif mod(iteration, roundtemp) == 0
        fprintf(repmat('\b', 1, barlength+3));
        fprintf(repmat('=', 1, floortemp));
        fprintf('>');
        fprintf(repmat(' ', 1, barlength-1-floortemp));
        fprintf(']\n\n'); 
    end


end
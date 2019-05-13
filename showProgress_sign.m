function showProgress_sign(iteration, maxvalue, division, barlength, progress_sign, progress_bar)
% Function that will fprintf a progression bar for any iterative algorithm.
%   -iteration: current iteration in the algorithm,
%   -maxvalue: maximum iteration in the algorithm,
%   -division (default 50): number of time progress bar will refresh,
%    division should be lower than maxvalue,
%   -barlength (default 50): length of progress bar,
%   -progress_sign (default '>'): displayed character for progress,
%   -progress_line (default '='): displayed line for progress.


    %% Default values:
    
    % Automatic values:
    if nargin == 5
        progress_bar = '=';
    elseif nargin == 4
        progress_bar = '=';
        progress_sign = '>';
    elseif nargin == 3
        progress_bar = '=';
        progress_sign = '>';
        barlength = 50;
    elseif nargin == 2
        progress_bar = '=';
        progress_sign = '>';
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
    % Check project_sign and progress_line are characters:
    if ~ischar(progress_sign) 
        error('progress_sign should be a character.')
    end
    progress_sign = char(progress_sign);
    lsign = length(progress_sign);
    if ~ischar(progress_bar)
        error('progress_line should be a character.')
    end
    progress_bar = char(progress_bar);
    if length(progress_bar) ~= 1
        error('progress_line should be of length 1.')
    end

    
    %% Initial empty progress bar:
    
    if iteration == 1
        fprintf('[');
        fprintf(progress_sign);
        fprintf(repmat(' ', 1, barlength-lsign));
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
        fprintf(repmat(progress_bar, 1, barlength));
        fprintf(']\n');
    elseif mod(iteration, roundtemp) == 0
        fprintf(repmat('\b', 1, barlength+3));
        fprintf(repmat(progress_bar, 1, floortemp));
        fprintf(progress_sign);
        fprintf(repmat('\b', 1, floortemp+lsign-barlength));
        fprintf(repmat(' ', 1, barlength-lsign-floortemp));
        fprintf(']\n\n'); 
    end


end
function [avg, avg_t] = average_signal(signal, times, period, style)
% This function averages signal, on a time frame equal to period.
% What it does is it interpolates signal from times to a multiple of the
% period, before averaging it over a period.
%
% Inputs:
% - signal [1D or 2D array]: signal to average. If several signals are
%   provided to be averaged, time must be the second axis.
% - times [1D array]: time at which each point from signal is sampled.
% - period [double]: time on which to average signal.
%
% Outputs:
% - avg [1D or 2D array]: averaged signal. If several signals were
%   provided, time is on the second axis.
% - avg_t [1D array]: times associated to averaged array.


    %% Rearange signal if necessary
    
    % If signal is a vector, reshape it is necessary
    if size(signal, 1) == 1 || size(signal, 2) == 1
        signal = reshape(signal, 1, length(signal));
    end
    
    % Get shape of signal
    [m, n] = size(signal);


    %% Find new time points and interpolate
    
    % Rearange times
    times = reshape(times, 1, n) - times(1); % first value is 0
    
    % Define new time vector
    period_closest = find(period <= times, 1); % find point closest to period
    % We are going to scale this index so that it gets to period
    dt = period / (period_closest-1); % find the new time increment after scaling
    t = times(1):dt:times(end);
    
    % Interpolate based on these points
    s = interp1(times, signal', t, 'spline')';
    if m == 1; s = reshape(s, 1, length(s)); end % otherwise it does not work with 1D arrays
    
    
    %% Average on a period
    
    % Define output
    avg = zeros(m, period_closest-1);
    
    % Fill output by looping for each period
    for p = 1:(period_closest-1)
        % Index of signal points
        ind = p:(period_closest-1):n;
        if nargin == 4 && isequal(style, 'median')
            avg(:, p) = median(s(:, ind), 2); 
        else
            avg(:, p) = sum(s(:, ind), 2) ./ length(ind);
        end
    end
    
    % Add first row as last row for periodicity
    avg = cat(2, avg, avg(:, 1));
    
    % Time array associated
    avg_t = t(1:period_closest);


end
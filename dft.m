function [f, A, phi, low_pass] = dft(signal, fs, varargin)
% This function performs a discrete Fourier transform of the given signal, 
% by projecting it onto a base of complex exponentials.
% NB: see https://www.youtube.com/watch?v=spUNpyF58BY which is an amazing
% representation of how Fourier series work.
%
% Inputs:
% - signal [1D or 2D array]: observed signal. If multiple signals are
%   provided, the time component should be along the 2nd axis.
% - fs [double]: sampling frequency (same for all). Note that the sampling 
%   theorem tells us that the sampling frequency fs must be superior or 
%   equal to 2 times the maximum frequency present in the signal (the 
%   Nyquist rate).
% - fcut [double, optional]: frequency for low pass filter. If none
%   provided, set to Nyquist frequency. 
% - fanalysis [double, optional]: if this frequency is provided, the
%   function is just going to project signal on 1 complex exponential with
%   frequency fanalysis, and return outputs based on that.
%
% Outputs:
% - f [1D or 2D array]: array of frequencies analyzed.
% - A [1D or 2D array]: array of amplitudes associated to frequencies.
% - phi [1D or 2D array]: array of phases associated to frequencies.
% - low_pass [function]: signal(s) rebuilt from Fourier coefficients, 
%   cutting frequency at fcut. Argument of function is time, and the index
%   of the subsignal in signal can also be provided.


    %% Check inputs
    
    % Input parser
    p = inputParser;
    addRequired(p, 'signal');
    addRequired(p, 'fs');
    addOptional(p, 'fcut', nan);
    addOptional(p, 'fanalysis', nan);
    parse(p, signal, fs, varargin{:});


    %% Parameters
    
    % Input signal
    if size(p.Results.signal, 1) == 1 || size(p.Results.signal, 2) == 1
        N = length(p.Results.signal);
        p.Results.signal = reshape(p.Results.signal, 1, N);
    else
        N = size(p.Results.signal, 2);
    end
    
    % DOES NOT RETURN THE SAME VALUE AS GENERAL ANALYSIS
    % If fanalysis is provided, work is easier
    if ~isnan(p.Results.fanalysis)
        % Output frequency
        f = p.Results.fanalysis;
        % Fourrier coefficient
        Xk = p.Results.signal * exp(-2*pi*1i*f .* (0:(N-1))');
        % Amplitude
        A = 2 * abs(Xk) / N;
        % Phase
        phi = angle(Xk);
        % Low pass
        if isnan(p.Results.fcut)
            fcut = f(end) + 0.001;
        else
            fcut = p.Results.fcut;
        end
        cut = (f < fcut);
        low_pass = @(varargin) low_pass_filter(A, f, phi, cut, varargin);
        % Exit function
        return
    end
    
    % Build frequencies output
    f = p.Results.fs / N .* (0:floor(N/2));
    
    
    %% Compute discrete Fourier transform
    
    % Build set of complex exponentials
    k = (0:(N-1)) .* ones(N); % index for Fourier coefficient
    n = (0:(N-1))' .* ones(N); % index for point considered
    expi = exp(-2*pi*1i/N .* k .* n);
    
    % Find Fourier coefficients
    Xk = p.Results.signal * expi;
    
    
    %% Get amplitudes and phases for these coefficients
    
    % Get module for complex number
    Amod = abs(Xk);
    % Find amplitude
    A = 2 * Amod(:, 1:floor(N/2)+1) / N;
    A(:, 1) = A(:, 1) / 2; % first amplitude is not x2
    
    % Get angle for complex number
    phimod = angle(Xk);
    phi = phimod(:, 1:floor(N/2)+1);
    
    
    %% Build low pass filter with cosines
    
    % Define cut frequency if not provided
    if isnan(p.Results.fcut)
        fcut = f(end) + 0.001;
    else
        fcut = p.Results.fcut;
    end
    
    % Build low pass filter based on fcut
    cut = (f < fcut);
    low_pass = @(varargin) low_pass_filter(A, f, phi, cut, varargin);
        
        
end

function out = low_pass_filter(A, f, phi, cut, varargin)
% This function is going to be low_pass filled with the right values for A,
% f, phi and cut.
% varargin is the time, and there is an optional parameter corresponding to
% the subsignal of the signal to low pass, so that all signal do not have
% to be computed at all times.

    if length(varargin{1}) == 1
        out = squeeze(sum(A(:, cut) .* cos(2*pi*f(:, cut).*reshape(varargin{1}{1}, 1, 1, length(varargin{1}{1})) + phi(:, cut)), 2));
    elseif length(varargin{1}) == 2
        out = squeeze(sum(A(varargin{1}{2}, cut) .* cos(2*pi*f(:, cut).*reshape(varargin{1}{1}, 1, 1, length(varargin{1}{1})) + phi(varargin{1}{2}, cut)), 2));
    end
    
end
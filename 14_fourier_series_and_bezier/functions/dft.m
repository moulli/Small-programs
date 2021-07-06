function [f, A, phi, low_pass] = dft(signal, fs, fcut)
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
%
% Outputs:
% - f [1D or 2D array]: array of frequencies analyzed.
% - A [1D or 2D array]: array of amplitudes associated to frequencies.
% - phi [1D or 2D array]: array of phases associated to frequencies.
% - low_pass [function]: signal(s) rebuilt from Fourier coefficients, 
%   cutting frequency at fcut. Argument of function is time.


    %% Parameters
    
    % Input signal
    if size(signal, 1) == 1 || size(signal, 2) == 1
        N = length(signal);
        signal = reshape(signal, 1, N);
    else
        N = size(signal, 2);
    end
    
    % Build frequencies output
    f = fs / N .* (0:floor(N/2));
    
    
    %% Compute discrete Fourier transform
    
    % Build set of complex exponentials
    k = (0:(N-1)) .* ones(N); % index for Fourier coefficient
    n = (0:(N-1))' .* ones(N); % index for point considered
    expi = exp(-2*pi*1i/N .* k .* n);
    
    % Find Fourier coefficients
    Xk = signal * expi;
    
    
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
    if nargin == 2
        fcut = f(end) + 0.001;
    end
    
    % Build low pass filter based on fcut
    cut = (f < fcut);
    low_pass = @(t) squeeze(sum(A(:, cut) .* cos(2*pi*f(:, cut).*reshape(t, 1, 1, length(t)) + phi(:, cut)), 2));
        
        
end
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
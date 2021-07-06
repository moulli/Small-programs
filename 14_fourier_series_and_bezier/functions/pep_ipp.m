function integrated = pep_ipp(a, b, q, k, alpha, beta)
% PEP_IPP stands for polynomial/exponential product integration par partie.
% This function integrate alpha*(t-k)^q * exp(beta*t) between a and b,
% using and analytical formula that I computed.
%
% Inputs:
% - a [double]: lower limit for integration.
% - b [double]: higher limit for integration.
% - q [integer]: polynomial power of t.
% - k [double]: subtracted to t before the power
% - alpha [complex]: coefficient before (t-k)^q.
% - beta [complex]: coefficient before t in exponential.
%
% Outputs:
% - integrated [complex]: result of Integration Par Partie.


    %% Compute integrated value
    
    % Define sum variable
    m = 1:(q+1);
    
    % Define intermediate matrix
    inm = (-1).^(m+1) .* factorial(q) ./ factorial(q+1-m) .* alpha ./ beta.^m ...
           .* ((b-k).^(q+1-m) .* exp(beta*b) - (a-k).^(q+1-m) .* exp(beta*a));
       
    % Sum over m to have result
    integrated = sum(inm);   
    
    % If beta is 0, computation is different
    if beta == 0
        integrated = alpha / (q+1) * ((b-k)^(q+1) - (a-k)^(q+1));
    end
    
    
end
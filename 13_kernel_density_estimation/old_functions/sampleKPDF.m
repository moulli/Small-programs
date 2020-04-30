function randkpdf = sampleKPDF(kpdf, n)

%% sampleKPDF returns points sampled from kpdf distribution
%
%  The kernel density estimation is a probability function, continuous over
%  the space specified. It allows to have the probability that a specified
%  point belongs to the space defined by coordinates. In order to get this
%  probability distribution, we use a simple centered standard gaussian
%  kernel. The only parameter is h, the bandwidth, which is a 'smoothening'
%  parameters: the highest it is, the more space kde will cover. The lowest
%  it is the more discrete the final distribution will be.
%  See: https://en.wikipedia.org/wiki/Kernel_density_estimation
%
%  Inputs:
%  - coordinates [(n x d) matrix]: n being the number of points to
%    consider, and d being the number of dimensions of the space.
%  - h [single/double number]: bandwidth (see above).
%
%  Outputs:
%  - kpdf [function]: kernel density estimation probability distribution.


    %% Check inputs
    
    % size of coordinates
    [~, d] = size(coordinates);
    
    % check bandwidth is a positive number
    if ~isequal(size(h), [1, 1]) || h <= 0
        error('bandwidth must be a positive number')
    end
    
    
    %% Build kernel density probability function
    
    % function for one point (kernel definition)
    kernel = @(x) mean( mvnpdf( (x-coordinates) / bandwidth ) / bandwidth.^d );    
    
    % function for many points
    randkpdf = @(x) cellfun( kernel, mat2cell( x, ones(size(x, 1), 1) ) );    
    
    
end
function out = SGfilterArray(data, rank, nleft, nright)
%% Function that performs Savitzky-Golay filtering on data.
%  The idea is inspired by the moving window averaging filtering algorithm,
%  which consists in taking, for each point of the data we want to
%  smoothen, the average value between it and its neighbours. Although in
%  this particular algorithm, we use a polynomial regression to do that.
%  Higher degree polynoms will have a better rendering for crests, and
%  therefore be more useful when it comes to analysing peaks.
%
%  PARAMETERS:
%  -'data': array of noisy values.
%  -'rank': rank of the polynom we want to fit.
%  -'nleft': number of points we want to take from the left.
%  -'nright'[optional]: same on the right. 
%  -'out': returns an array containing the filtered values.



    %% Parameters and errors:
    
    % Definition of window:
    if nargin == 3
        nright = nleft;
    end
    window = nleft + nright + 1;  
    % Reorganising data:
    [idata, jdata] = size(data);
    if idata ~= 1 && jdata ~= 1
        error('Please provide data under the form of an array')
    elseif length(data) < window
        error('Please provide a window shorter than number of data')
    end
    ldata = length(data);
    data = reshape(data, ldata, 1);
    % Definition of out:
    out = zeros(1, ldata);
    
    
    
    %% Filtering:
    
    for i = 1:ldata
        % Finding window for point of interest:
        nlefttemp = min([i-1, nleft]) + (nright-ldata+i) * (ldata-i < nright);
        nrighttemp = min([ldata-i, nright]) + (nleft-i+1) * (i-1 < nleft);
        % Taking points in the window:
        winpt = data(i-nlefttemp:i+nrighttemp);
        % Moore-Penrose for polynomial linear regression:
        absc = (-nlefttemp:nrighttemp)';
        degree = 0:rank;
        Xpol = absc .^ degree;
        weights = (Xpol' * Xpol) \ Xpol' * winpt;
        % Inserting new point in array:
        out(i) = Xpol(nlefttemp+1, :) * weights;
    end
    


end
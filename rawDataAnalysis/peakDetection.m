function out = peakDetection(data, approach, param, space)
%% Function that detects the maximum peaks in an array of values.
%  The algorithm computes all the negative derivatives and deduces which
%  values constitute a local maximum. Based on parameters given by the
%  user, the function then returns all the peaks in the data. Filtering the
%  data, for example with fSavitzkygolayFilter, is recommended.
%
%  PARAMETERS:
%  -'data': array of values to get peaks from.
%  -'approach': string in {'abs', 'rel', 'pro', 'num'}; it determines the 
%   approach to keep peaks.
%  -'param': either the value starting which we keep the peaks (absolute
%   value or relative to max(data), or the proportion of peaks we keep.
%  -'space'[optional]: index between two points at which we only consider
%   one peak, chosen to be the maximum between the n peaks involved. We
%   also get rid of the small maxima between two minima with this
%   technique.



    %% Parameters and errors:
    
    % Definition of space parameter:
    if nargin == 3
        space = 1;
    end 
    % Reorganising data:
    [idata, jdata] = size(data);
    if idata ~= 1 && jdata ~= 1
        error('Please provide data under the form of an array')
    end
    ldata = length(data);
    data = reshape(data, ldata, 1);
    % Approach to keep peaks:
    if sum(approach == ['abs'; 'rel'; 'pro'; 'num']) == 0
        error('Approach must be "abs", "rel", "pro", or "num"')
    end
    % Definition of out:
    out = zeros(1, ldata);
    
    
    
    %% Finding all local maxima & minima in the data:
    
    % Positive and negative differences:
    diffv = data(2:end)-data(1:end-1);
    diffv = (diffv > 0) - (diffv < 0);
    % Deducing maxima and minima:
    evolution = [0; (diffv(1:end-1)-diffv(2:end))/2; 0];
    maxe = find(evolution == 1);
    mine = find(evolution == -1);
    
    
    
    %% Using space parameter to get rid of some peaks:
    
    
    
    %% Defining peaks to keep based on approach:
    if approach == 'abs'
        out = maxe(data(maxe) >= param);
    elseif approach == 'rel'
        out = maxe(data(maxe) >= max(data) * param);
    elseif approach == 'pro'
        lmax = length(maxe);
        lkeep = ceil(param * lmax);
        [~, indmaxe] = sort(data(maxe), 'descend');
        out = maxe(indmaxe(1:lkeep));
    else
        [~, indmaxe] = sort(data(maxe), 'descend');
        lkeep = min(param, length(indmaxe));
        out = maxe(indmaxe(1:lkeep));
    end    
    


end
function [Vrid, varargout] = static_ridNeurons(Vin, ridvalues, varargin)

%% Function that will get rid of values outside provided interval, and NaNs.
%
%  This function will delete values in Vin who are higher than the first 
%  argument of ridvalues, and lower than the second argument of ridvalues. 
%  Ultimately the same ROWS deletion will be applied to any input provided 
%  additionally to these two mandatory inputs.
%
%
%% Inputs:
%
%  --Vin: array of values we want to analyse.
%  --ridvalues: array composed of 2 values. First one is lower bound of
%    deletion interval and second one is upper bound of deletion interval. 
%    Note that if ridvalues(1) > ridvalues(2), algorithm will delete values
%    that are outside interval.
%  --varargin: arrays or matrices with as many rows as there are values in
%    Vin. Rows corresponding to deteted values from Vin will also be
%    deleted in varargin.
%
%
%% Outputs:
%
%  --Vrid: Vin without deleted values.
%  --varargout: varargin without deleted row values. Contains indexes of
%    removed values as last cell.



    %% Initialization:
    
    % Checking Vin is a vector:
    if ~isvector(Vin)
        error('Please provide Vin as a vector.')
    end
    % Checking ridvalues is a vector of two values:
    if ~isvector(ridvalues) || (isvector(ridvalues) && length(ridvalues) ~= 2)
        error('Please provide ridvalues as a two values vector.')
    else
        if ridvalues(1) < ridvalues(2)
            inside = 1;
        elseif ridvalues(1) > ridvalues(2)
            inside = 0;
        else
            error('Please provide ridvalues with different values.')
        end
    end
    % Building varargout based on varargin:
    lvar = length(varargin);
    varargout = cell(length(varargin), 1);
    
    
    
    %% Finding values to delete:
    
    % Getting values from ridvalues:
    if inside == 1
        zerosN = (ridvalues(1) <= Vin & Vin <= ridvalues(2));
    elseif inside == 0
        zerosN = (Vin <= ridvalues(2) | ridvalues(1) <= Vin);
    end
    % Adding nan values:
    vinan = isnan(Vin);
    zerosN = (zerosN | vinan);
    
    
    
    %% Now deleting these values from Vin and varargin:
    
    % From Vin to Vrid:
    Vrid = Vin;
    Vrid(zerosN) = [];
    % From varargin to varargout:
    for i = 1:lvar
        vartemp = varargin{i};
        vartemp(zerosN, :) = [];
        varargout{i} = vartemp;
    end    
    % Adding removed values:
    varargout{lvar+1} = find(zerosN);


end
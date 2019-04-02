function Ccolor = static_corr2col(Ccorrelation, varargin)

%% Function that adapts higher correlation to darker color.
%
%
%% Input:
%
%  --Ccorrelation: array of correlation, ranging from -1 to 1.
%  --varargin: can delete autoscale, by writing 'autoscale', 'off'. By
%    default, all values are autoscaled, meaning whatever the input, the 
%    output will range from 0 to 1. 
%
%
%% Output:
%
%  --Ccolor: color matrix, with same number of rows as Ccorrelation length,
%    and 3 columns.



    %% Initialization:
    
    % Confirm Ccorrelation is a vector:
    if ~isvector(Ccorrelation)
        error('Please provide Ccorrelation as a vector.')
    end
    % Get autoscale information:
    autoscale = 1;
    if nargin ~= 1
        if nargin ~= 3
            error('Please provide right attribute.')
        elseif string(varargin{1}) ~= 'autoscale'
            error('Input is not an attribute of the function.')
        elseif sum(string(varargin{2}) == ["true", "false"]) == 0
            error('Input does not correspond to attribute setting.')
        elseif string(varargin{2}) == "false"
            autoscale = 0;
        end
    end
    
    
    
    %% Defining Ccolor matrix:
   
    % Inverting values:
    Ccolorn1 = 1 - Ccorrelation.*(Ccorrelation > 0);
    Ccolorn2 = 1 + Ccorrelation.*(Ccorrelation < 0);
    % Autoscaling if required:
    if autoscale 
        if length(unique(Ccolorn1)) ~= 1
            Ccolorn1 = (Ccolorn1 - min(Ccolorn1)) / (max(Ccolorn1) - min(Ccolorn1));
        end
        if length(unique(Ccolorn2)) ~= 1
            Ccolorn2 = (Ccolorn2 - min(Ccolorn2)) / (max(Ccolorn2) - min(Ccolorn2));
        end
    end
    Ccolor = [Ccolorn2, Ccolorn1.*(Ccorrelation > 0)+Ccolorn2.*(Ccorrelation < 0), Ccolorn1];


end
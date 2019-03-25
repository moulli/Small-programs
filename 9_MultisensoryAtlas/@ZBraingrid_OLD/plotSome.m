function plotSome(obj, subset, varargin)

%% Function that plots averaged correlation of a subset in the ZBraingrid class.
%
%  Takes in input the properties of the object, and plots a 3D
%  representation of the brain, with the correlations averaged over a
%  subset of dataset for all the points in the grid.
% 
%
%% Inputs:
%
%  --obj: references the object this methods is attached to.
%  --subset: way of selecting the subset to be plot. Can either be a vector
%    of the datasets the user wants to plot, or a string, which will
%    trigger a research in the comment property of the ZBraingrid object.
%  --varargin: three optional inputs can be added as inputs to the plotAll
%    method. First 'ridvalues', associated to a vector of two values gets
%    rid of some neurons (cf. aux_ridNeurons). Second 'autoscale', allows
%    to scale colors in plot (cf. aux_corr2col). Third 'MarkerSize'
%    (default 30) sets the size of the points in the scatter plot.
%    'intercept' allows to only plot points that are shared between
%    datasets, and 'autocrop' on 'off' allows to adapt color scale after
%    ridvalues.



    %% Initialization:
    
    % Checking inputs:
    intext = ["ridvalues", "autoscale", "MarkerSize", "intercept", "autocrop"];
    ridindication = 0;
    scaleindication = 0;
    markersize = 30;
    intindication = 0;
    cropindication = 0;
    if nargin > 2
        if sum(nargin == [4, 6, 8, 10]) == 0
            error('Input arguments must come by pairs.')
        else
            for i = 1:((nargin-2)/2)
                optindex = find(varargin{(i-1)*2+1} == intext);
                switch optindex
                    case 1
                        ridindication = 1;
                        ridvalues = varargin{2*i};
                    case 2
                        scaleindication = 1;
                        autoscale = varargin{2*i};
                    case 3
                        markersize = varargin{2*i};
                    case 4
                        if lower(varargin{2*i}) == true
                            intindication = 1;
                        elseif lower(varargin{2*i}) ~= false
                            error('Attribute associated to intercept is wrong.')
                        end
                    case 5
                        if varargin{2*i} == false
                            cropindication = 1;
                        elseif lower(varargin{2*i}) ~= true
                            error('Attribute associated to CropScale is wrong.')
                        end
                    otherwise
                        error('Input arguments must be ridvalues, autoscale, or MarkerSize.')
                end
            end
        end
    end
    % Chosing subsets:
    if ~isvector(subset)
        error('Please provide subset as a vector or as a string.')
    end
    if isnumeric(subset)
        checksub = sum((1:length(obj.names))' == reshape(subset, 1, length(subset)));
        if ~isequal(checksub, ones(size(checksub)))
            error('Please provide subsets in range of available sets.')
        else
            keepsub = subset;
        end
    else
        keepsub = [];
        for i = 1:length(obj.comments)
            if ~isempty(strfind(obj.comments{i}, subset))
                keepsub = [keepsub, i];
            end
        end
    end
    if isempty(keepsub)
        error('Subset to plot is empty.')
    end
    
    
    
    %% Creating meshgrid:
    
    xmesh = (obj.xgrid(1:end-1)+obj.xgrid(2:end)) / 2;
    ymesh = (obj.ygrid(1:end-1)+obj.ygrid(2:end)) / 2;
    zmesh = (obj.zgrid(1:end-1)+obj.zgrid(2:end)) / 2;
    [Xmesh, Ymesh, Zmesh] = meshgrid(xmesh, ymesh, zmesh);
    Xmesh = permute(Xmesh, [2, 1, 3]);
    Ymesh = permute(Ymesh, [2, 1, 3]);
    Zmesh = permute(Zmesh, [2, 1, 3]);
    % Grid coordinates:
    g_coord = [Xmesh(:), Ymesh(:), Zmesh(:)];
    % Getting rid of unnecessary neurons:
    Cgridm = mean(obj.Zcorrelations(:, :, :, keepsub), 4);
    Cgridm = Cgridm(:);
    % If intercept, just keeping neurons in all datasets:
    if intindication == 1
        Itemp = all(obj.Zneuron_number, 4);
        Cgridm = Cgridm(Itemp(:));
        g_coord = g_coord(Itemp(:), :);
    end
    % Defining color scaling:
    if scaleindication == 1
        Ccolor = ZBraingrid.static_corr2col(Cgridm, 'autoscale', autoscale);
    else
        Ccolor = ZBraingrid.static_corr2col(Cgridm);
    end
    % Getting rid of values:
    if ridindication == 1
        [Cgridmlayed, g_coord, Ccolor] = ZBraingrid.static_ridNeurons(Cgridm, ridvalues, g_coord, Ccolor);
    else
        Cgridmlayed = Cgridm;
        zerocor = (Cgridmlayed == 0);
        Cgridmlayed(zerocor) = [];
        g_coord(zerocor, :) = [];
        Ccolor(zerocor, :) = [];
    end
    
    
    
    %% Plotting:
    
    % Defining color scaling if no crop:
    if cropindication == 1
        if scaleindication == 1
            Ccolor = ZBraingrid.static_corr2col(Cgridmlayed, 'autoscale', autoscale);
        else
            Ccolor = ZBraingrid.static_corr2col(Cgridmlayed);
        end
    end
    
    % Plotting:
    figure
    scatter3(g_coord(:, 1), g_coord(:, 2), g_coord(:, 3), markersize, Ccolor, 'filled')
    axis equal
    title('Grid inside Zbrain with correlation', 'Interpreter', 'latex')
    xlabel('x-axis', 'Interpreter', 'latex')
    ylabel('y-axis', 'Interpreter', 'latex')
    zlabel('z-axis', 'Interpreter', 'latex')
           

end
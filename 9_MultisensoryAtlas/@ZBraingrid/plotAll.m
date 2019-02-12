function plotAll(obj, varargin)

%% Function that plots averaged correlation in the ZBraingrid object.
%
%  Takes in input the properties of the object, and plots a 3D
%  representation of the brain, with the correlations averaged over every
%  dataset for all the points in the grid.
% 
%
%% Inputs:
%
%  --obj: references the object this methods is attached to.
%  --varargin: three optional inputs can be added as inputs to the plotAll
%    method. First 'ridvalues', associated to a vector of two values gets
%    rid of some neurons (cf. aux_ridNeurons). Second 'autoscale', allows
%    to scale colors in plot (cf. aux_corr2col). Third 'MarkerSize'
%    (default 30) sets the size of the points in the scatter plot.



    %% Initialization:
    
    % Indication:
    fprintf('\n\nLaunching function plotAll, attribute of ZBraingrid class.\n');
    % Checking inputs:
    intext = ["ridvalues", "autoscale", "MarkerSize"];
    ridindication = 0;
    scaleindication = 0;
    markersize = 30;
    if nargin > 1
        if sum(nargin == [3, 5, 7]) == 0
            error('Input arguments must come by pairs.')
        else
            for i = 1:((nargin-1)/2)
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
                    otherwise
                        error('Input arguments must be ridvalues, autoscale, or MarkerSize.')
                end
            end
        end
    end
    
        
%     if nargin > 1
%         switch nargin
%             case 3
%                 switch lower(varargin{1})
%                     case intext(1)
%                         ridindication = 1;
%                         ridvalues = varargin{2};
%                     case intext(2)
%                         scaleindication = 1;
%                         autoscale = varargin{2};
%                     otherwise
%                         error('Input arguments must be ridvalues or autoscale.')
%                 end
%             case 5
%                 switch lower(varargin{1})
%                     case intext(1)
%                         ridindication = 1;
%                         ridvalues = varargin{2};
%                         switch lower(varargin{3})
%                             case intext(2)
%                                 scaleindication = 1;
%                                 autoscale = varargin{4};
%                             otherwise
%                                 error('Input arguments must be ridvalues or autoscale.')
%                         end
%                     case intext(2)
%                         scaleindication = 1;
%                         autoscale = varargin{2};
%                         switch lower(varargin{3})
%                             case intext(1)
%                                 ridindication = 1;
%                                 ridvalues = varargin{4};
%                             otherwise
%                                 error('Input arguments must be ridvalues or autoscale.')
%                         end
%                     otherwise
%                         error('Input arguments must be ridvalues or autoscale.')
%                 end
%         end
%     end
    
    
    
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
    Cgridm = mean(obj.Zcorrelations, 4);
    if ridindication == 1
        [Cgridmlayed, g_coord] = ZBraingrid.aux_ridNeurons(Cgridm(:), ridvalues, g_coord);
    else
        Cgridmlayed = Cgridm(:);
    end
    
    
    
    %% Plotting:
    
    % Defining color scaling:
    if scaleindication == 1
        Ccolor = ZBraingrid.aux_corr2col(Cgridmlayed, 'autoscale', autoscale);
    else
        Ccolor = ZBraingrid.aux_corr2col(Cgridmlayed);
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
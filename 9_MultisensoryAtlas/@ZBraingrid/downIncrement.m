function onew = downIncrement(obj, new_increment)

%% Function that created a new object with given increment in the ZBraingrid class.
%
%  Takes in input the ZBraingrid object, checks new_increment is lower to
%  current increment, and computes properties for new object.
% 
%
%% Input:
%
%  --obj: references the object this methods is attached to.
%  --new_increment: increment for new object. Must be inferior to former 
%    increment.
%
%
%% Output:
%
%  --onew: new ZBraingrid object. 



    %% Initialization:
    
    % Checking number of arguments:
    if nargin ~= 2
        error('Please provide ZBraingrid object and new increment.')
    end
    % Manipulating increment:
    if ~isvector(new_increment)
        error('Please provide new increment as a vector.')
    elseif length(new_increment) == 1
        new_increment = new_increment .* ones(1, 3);
    end
    if length(new_increment) ~= 3
        error('Please provide new increment as a scalar or as a 1x3 vector.')
    end
    % Checking new increment:
    if all(new_increment <= obj.increment)
        error('Please provide increment lower than current increment.')
    else
        onew = ZBraingrid(obj.method, new_increment);
    end
    
    
    
    %% Filling same properties:
    
    onew.names = obj.names;
    onew.paths = obj.paths;
    onew.comments = obj.comments;
    onew.Zcorvect = obj.Zcorvect;
    
    
    
    %% Filling computed properties:
    
    [lx, ly, lz] = size(obj.Zneurons);
    for ix = 1:lx
        for iy = 1:ly
            for iz = 1:lz
%                 % Getting neurons in this part of grid:
%                 xtemp = (onew.xgrid(ix) <= coord_in(:, 1) & coord_in(:, 1) < onew.xgrid(ix+1));
%                 ytemp = (onew.ygrid(iy) <= coord_in(:, 2) & coord_in(:, 2) < onew.ygrid(iy+1));
%                 ztemp = (onew.zgrid(iz) <= coord_in(:, 3) & coord_in(:, 3) < onew.zgrid(iz+1));
%                 Ttemp = find(xtemp & ytemp & ztemp);
%                 % Now updating zbraingrid object:
%                 Zneurons_temp{ix, iy, iz} = Ttemp;
%                 cortemp = mean(cor_in(Ttemp));
%                 if isnan(cortemp); cortemp = 0; end
%                 Zcorrelations_temp(ix, iy, iz) = mean(cortemp);
%                 Zneuron_number_temp(ix, iy, iz) = length(Ttemp);
            end
        end
    end


end
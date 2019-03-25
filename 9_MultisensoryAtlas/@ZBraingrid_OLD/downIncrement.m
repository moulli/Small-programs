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
    
    % Getting necessary information:
    [lx, ly, lz, ~] = size(onew.Zneurons);
    ld = size(obj.Zneurons, 4);
    xinc = (obj.xgrid(1:end-1) + obj.xgrid(2:end)) ./ 2;
    yinc = (obj.ygrid(1:end-1) + obj.ygrid(2:end)) ./ 2;
    zinc = (obj.zgrid(1:end-1) + obj.zgrid(2:end)) ./ 2;
    % Concatenating matrices:
    for id = 1:(ld-1)
        onew.Zneurons = cat(4, onew.Zneurons, cell(lx, ly, lz, 1));
        onew.Zcorrelations = cat(4, onew.Zcorrelations, zeros(lx, ly, lz, 1));
        onew.Zneuron_number = cat(4, onew.Zneuron_number, zeros(lx, ly, lz, 1));
    end
    % New neuron assignment:
    for ix = 1:lx
        for iy = 1:ly
            for iz = 1:lz
                for id = 1:ld
                    [ix, iy, iz, id];
                    % Getting neurons in this part of grid:
                    xtemp = (onew.xgrid(ix) <= xinc & xinc < onew.xgrid(ix+1));
                    ytemp = (onew.ygrid(iy) <= yinc & yinc < onew.ygrid(iy+1));
                    ztemp = (onew.zgrid(iz) <= zinc & zinc < onew.zgrid(iz+1));
                    Ttemp = obj.Zneurons(xtemp, ytemp, ztemp, id);
                    onew.Zneurons{ix, iy, iz, id} = sort(cell2mat(Ttemp(:)));
                    % Computing rest of object:
                    cortemp = mean(onew.Zcorvect{id}(onew.Zneurons{ix, iy, iz, id}));
                    if isnan(cortemp); cortemp = 0; end
                    onew.Zcorrelations(ix, iy, iz, id) = cortemp;
                    onew.Zneuron_number(ix, iy, iz, id) = length(Ttemp);
                end
            end
        end
    end
    onew.gridsize = size(onew.Zcorrelations);


end
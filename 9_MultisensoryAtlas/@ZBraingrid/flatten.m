function onew = flatten(obj, opt_comment)

%% Function that flattens object along 4th dimension in the ZBraingrid class.
%
%  Takes in input the ZBraingrid object, checks information corresponds
%  along properties, and averages data. Zneurons is set to a cell of 0,
%  since it has no sense to keep this property. Zcorrelations is averaged,
%  and Zneuron_number is summed, all these processes along the 4th
%  dimension, that is across all datasets. names, paths and comments are
%  set to default values.
% 
%
%% Input:
%
%  --obj: references the object this methods is attached to.
%  --opt_comment: optional comment for new ZBraingrid object.
%
%
%% Output:
%
%  --onew: new ZBraingrid object, with a Zneurons property set to a cell of
%    0 to make sense. names, paths and comments are set to default values. 



    %% Initialization:
    
    % Creating new ZBraingrid object:
    onew = ZBraingrid(obj.method, obj.increment);
    % Adding information:
    onew.names = ["Flattened ZBraingrid object from " + string(length(obj.names)) + " datasets."];
    onew.paths = "No path for flattened ZBraingrid objects";
    if nargin == 2
        onew.comments = string(opt_comment);
    elseif length(unique(obj.comments)) == 1
        onew.comments = unique(obj.comments);
    else
        onew.comments = "No comment";
    end
    
    
    
    %% Computing Zcorrelations and Zneuron_number:
    
    onew.Zcorrelations = mean(obj.Zcorrelations, 4);
    onew.Zneuron_number = sum(obj.Zneuron_number, 4);
    onew.gridsize = size(onew.Zcorrelations);
    
    
    
    %% Deleting Zcorvect and Zneurons:
    
    onew.Zcorvect = {};
    onew.Zneurons = cell(size(onew.Zcorrelations));


end
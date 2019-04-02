function cleanDuplicates(obj)

%% Function that cleans duplicates (if they exist) in the ZBraingrid class.
%
%  Takes in input the ZBraingrid object, checks if there is a redundant
%  name, and if so, checks if the data associated to this duplicate are the
%  same. If so, the latest entry is deleted from the object. If not, user
%  is informed that a name exists several time, without concording data.
% 
%
%% Input:
%
%  --obj: references the object this methods is attached to.



    %% Finding duplicates:
    
    % Indication:
    fprintf('\nLaunching function cleanDuplicates, attribute of ZBraingrid class.\n');
    % Getting names:
    name = obj.names;
    % To keep vector:
    to_keep = ones(length(name), 1);
    % Checking if name is in double:
    same_mat = tril(string(name) == string(name)', -1);
    same_vect = sum(same_mat, 2);
    for i = length(same_vect):-1:1
        similar = find(same_mat(i, :) == 1);
        for j = 1:length(similar)
            simtemp = similar(j);
            if isequal(obj.paths(i), obj.paths(simtemp)) && isequal(obj.Zcorvect(i), obj.Zcorvect(simtemp)) && ...
                    isequal(obj.Zneurons(i), obj.Zneurons(simtemp)) && isequal(obj.Zcorrelations(i), obj.Zcorrelations(simtemp)) && ...
                    isequal(obj.Zneuron_number(i), obj.Zneuron_number(simtemp))
                to_keep(i) = 0;
                fprintf('Found a redundant dataset, between subsets %.0f and %.0f.\nDeleting subset %.0f.\n', [simtemp, i, i]);
                break
            else
                fprintf('Found a redundant name, between subsets %.0f and %.0f.\nNo subset deleted.\n', [simtemp, i]);
            end
        end
    end
    if isequal(to_keep, ones(length(name), 1))
        fprintf('No duplicate found.\n');
    end
    
    
    
    %% Removing duplicates:
    
    to_keep = (to_keep == 1);
    obj.names = obj.names(to_keep);
    obj.paths = obj.paths(to_keep);
    obj.comments = obj.comments(to_keep);
    obj.Zcorvect = obj.Zcorvect(to_keep);
    obj.Zneurons = obj.Zneurons(:, :, :, to_keep);
    obj.Zcorrelations = obj.Zcorrelations(:, :, :, to_keep);
    obj.Zneuron_number = obj.Zneuron_number(:, :, :, to_keep);
    obj.gridsize = size(obj.Zcorrelations);


end
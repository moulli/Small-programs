function [B, BC, P] = bezier_from_svg(path_to_file)
% This function extracts points and Bezier parameters from .svg file.
% !!! PATH MUST BE A COMPOSITION OF CUBIC BEZIER CURVES !!!
% It first import the file, and then extract the points responsible for the
% path. Using these points, it builds a set of parametric bezier curves,
% with each set corresponding to 4 points.
% See: https://www.youtube.com/watch?v=pnYccz1Ha34 and 
% https://css-tricks.com/svg-path-syntax-illustrated-guide/ for a full
% explaination.
% 
% Inputs:
% - path_to_file [string]: path to svg file.
%
% Outputs:
% - B [cell]: cell of parametric functions, one for each set of 4 points. 
% - BC [cell]: cell of coefficients. A cubic Bezier curve can be rewritten 
%   in the form a + b*t + c*t^2 + d*t^3. For each cell, there is the array 
%   [a, b, c, d].
% - P [cell]: cell returning the points in the svg path. Cell is organised 
%   so that the number of lines is the number of Bezier curves, and the 
%   number of columns is 4, i.e. the number of points to parametrize a 
%   cubic Bezier curve.


    %% Import svg file
    
    S = importdata(path_to_file);
    
    
    %% Extract points from svg
    
    % Find path location and therefore path
    for l = 1:length(S)
        if ~isempty(regexp(S{l}, 'path', 'once'))
            % Path line
            path = S{l+2};
            % Split
            path_split = split(path);
            % Parameters
            beginning = false;
            count = -1;
            relative = false;
            % Loop over split cell
            for sp = 1:length(path_split)
                % If first point
                if ~isempty(regexp(path_split{sp}, 'm', 'once'))
                    beginning = true;
                % If cubic bezier relative
                elseif ~isempty(regexp(path_split{sp}, 'c', 'once'))
                    relative = true;
                % If cubic bezier absolute
                elseif ~isempty(regexp(path_split{sp}, 'C', 'once'))
                    relative = false;
                % If end of path
                elseif ~isempty(regexp(path_split{sp}, 'z', 'once'))
                    % Remove last line and break
                    P = P(1:size(P, 1)-1, :);
                    break
                % Else add new point
                else %if ~isempty(path_split{sp}) % first one is empty
                    % Get point and break if not a number
                    point = path_split{sp};
                    point_split = split(point, ',');
                    try
                        point_num = [str2double(point_split{1}), str2double(point_split{2})];
                    catch
                        continue
                    end
                    % If first point start cell
                    if beginning
                        P = {};
                        P{1, 1} = point_num;
                        beginning = false;
                    % Else fill cell
                    else
                        % Find indexes
                        ind_bez = floor(count/3) + 1;
                        ind_pt = mod(count, 3) + 2;
                        % Add point in relative to former point
                        if relative
                            point_num = P{ind_bez, 1} + point_num;
                        end
                        % Add point to cell
                        P{ind_bez, ind_pt} = point_num;
                        if ind_pt == 4
                            P{ind_bez+1, 1} = point_num;
                        end
                    end
                    % Add count
                    count = count + 1;
                end
            end
        end
    end

    % Turn y-axis around (Inkscape and Matlab have a different y-axis orientation
    for i = 1:size(P, 1)
        for j = 1:size(P, 2)
            P{i, j}(2) = -P{i, j}(2);
        end
    end
    
    
    %% Define parametric output and coefficients
    
    % Parametric output
    B = cell(size(P, 1), 1);
    for p = 1:size(P, 1)
        B{p} = @(t) (1-t).^3.*P{p, 1} + 3*(1-t).^2.*t.*P{p, 2} + 3*(1-t).*t.^2.*P{p, 3} + t.^3.*P{p, 4};
    end
    
    % Bezier parametric coefficients
    BC = cell(size(P, 1), 1);
    for p = 1:size(P, 1)
        BC{p} = [P{p, 1};
                 3*(P{p, 2}-P{p, 1});
                 3*(P{p, 3}-2*P{p, 2}+P{p, 1});
                 P{p, 4}-3*P{p, 3}+3*P{p, 2}-P{p, 1}];
    end    
    
    
end
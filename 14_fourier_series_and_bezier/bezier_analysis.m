clear; close all; clc


%% Get SVG points

% Import file
file = 'bezier_svg_2.svg';
S = importdata(file);

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

% Turn y-axis around
for i = 1:size(P, 1)
    for j = 1:size(P, 2)
        P{i, j}(2) = -P{i, j}(2);
    end
end


% % Find path location and therefore path
% for l = 1:length(S)
%     if ~isempty(regexp(S{l}, 'path', 'once'))
%         % Path line
%         path = S{l+2};
%         % Split
%         path_split = split(path);
%         % Loop over split cell
%         for sp = 1:length(path_split)
%             % First point
%             if ~isempty(regexp(path_split{sp}, 'm', 'once'))
%                 % Define points cell
%                 P = {};
%                 % Get point
%                 point = path_split{sp+1};
%                 point_split = split(point, ',');
%                 P{1, 1} = [str2double(point_split{1}), str2double(point_split{2})];
%             % Other points for cubic bezier lines
%             elseif ~isempty(regexp(path_split{sp}, 'c', 'once'))
%                 % Get points
%                 for spc = (sp+1):length(path_split)
%                     % Break if end of segment
%                     if ~isempty(regexp(path_split{spc}, 'z', 'once')) || ~isempty(regexp(path_split{spc}, 'C', 'once'))
%                         break
%                     end
%                     point = path_split{spc};
%                     point_split = split(point, ',');
%                     ind_bez = floor((spc-sp-1)/3) + 1;
%                     ind_pt = mod(spc-sp-1, 3) + 2;
%                     % Add point in relative to former point
%                     Pnew = P{ind_bez, 1} + [str2double(point_split{1}), str2double(point_split{2})];
%                     % Add point to cell
%                     P{ind_bez, ind_pt} = Pnew;
%                     if ind_pt == 4
%                         P{ind_bez+1, 1} = Pnew;
%                     end
%                 end
%             elseif ~isempty(regexp(path_split{sp}, 'C', 'once'))
%                 % Get points
%                 for spc = (sp+1):length(path_split)
%                     % Break if end of segment
%                     if ~isempty(regexp(path_split{spc}, 'z', 'once')) || ~isempty(regexp(path_split{spc}, 'c', 'once'))
%                         break
%                     end
%                     point = path_split{spc};
%                     point_split = split(point, ',');
%                     ind_bez = floor((spc-sp-1)/3) + 1;
%                     ind_pt = mod(spc-sp-1, 3) + 2;
%                     % Add point in relative to former point
%                     Pnew = [str2double(point_split{1}), str2double(point_split{2})];
%                     % Add point to cell
%                     P{ind_bez, ind_pt} = Pnew;
%                     if ind_pt == 4
%                         P{ind_bez+1, 1} = Pnew;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % Remove last line 
% P = P(1:size(P, 1)-1, :);


%% Parametrize

% Parameter
t = (0:0.01:1)';
% Cubic bezier
cubic_bezier = @(t, p0, p1, p2, p3) (1-t).^3.*p0 + 3*(1-t).^2.*t.*p1 + 3*(1-t).*t.^2.*p2 + t.^3.*p3;
% Full path
Pt = zeros(0, 2);
for i = 1:size(P, 1)
    Pt = cat(1, Pt, cubic_bezier(t, P{i, 1}, P{i, 2}, P{i, 3}, P{i, 4}));
end

figure
subplot(5, 1, 1:3)
hold on
for i = 1:size(P, 1)
    for j = 1:size(P, 2)
        scatter(P{i, j}(1), P{i, j}(2))
    end
end
scatter(Pt(:, 1), Pt(:, 2), '.')
axis equal
subplot(5, 1, 4)
plot(Pt(:, 1), '.')
subplot(5, 1, 5)
plot(Pt(:, 2), '.')


%% Integrate polynoms

% Check it works with an approximation
a = 0;
b = 5;
q = 3;
k = 1.5;
alpha = 0.9;
beta = 3*pi*1i;
out = ipp(a, b, q, k, alpha, beta);
fprintf('Analytical result is %f + %fj \n', real(out), imag(out));
% Approximate with small surfaces
dx = 0.000001;
t = a:dx:b;
out_comp = sum(alpha.*(t-k).^q.*exp(beta.*t))*dx;
fprintf('Computational result is %f + %fj \n', real(out_comp), imag(out_comp));

% Define bezier coefficients
bezier_coeffs = @(p0, p1, p2, p3) [p0; 3*(p1-p0); 3*(p2-2*p1+p0); (p3-3*p2+3*p1-p0)];
bc = bezier_coeffs([0, 0], [0, 1], [1, 0], [1, 1]);
t = (0:0.001:1)';
temp = bc(1, :) + bc(2, :).*t + bc(3, :).*t.^2 + bc(4, :).*t.^3;
temp2 = cubic_bezier(t, [0, 0], [0, 1], [1, 0], [1, 1]); % ok it's the same formula

% Integrate
N = 25;
cn = zeros(2*N+1, 1);
T = size(P, 1);
for n = -N:N
    for p = 1:size(P, 1)
        bc = bezier_coeffs(P{p, 1}, P{p, 2}, P{p, 3}, P{p, 4});
        for bci = 1:size(bc, 1)
            a = p-1; b = p; q = bci-1;
            alpha = bc(bci, 1)+1i*bc(bci, 2);
            beta = -2*pi/T*n*1i;
            cn(n+N+1) = cn(n+N+1) + ipp(a, b, q, p-1, alpha, beta)/T;
        end
    end        
end

% Rebuild signal
t = (0:0.01:T)';
lent = length(t);
signal = exp(2*pi/T*(-N:N)*1i.*t) * cn;

figure
subplot(5, 1, 1:3)
hold on
scatter(Pt(:, 1), Pt(:, 2), 3, '.')
plot(signal)
axis equal
subplot(5, 1, 4)
hold on
plot(linspace(0, T, size(Pt, 1)), Pt(:, 1), '.')
plot(t, real(signal))
subplot(5, 1, 5)
hold on
plot(linspace(0, T, size(Pt, 1)), Pt(:, 2), '.')
plot(t, imag(signal))




        
                
function integrated = ipp(a, b, q, k, alpha, beta)
% This function integrate alpha*(t-k)^q * exp(beta*t) between a and b.
%
% Inputs:
% - a [double]: lower limit for integration.
% - b [double]: higher limit for integration.
% - q [integer]: polynomial power of t.
% - k [double]: subtracted to t before the power
% - alpha [complex]: coefficient before (t-k)^q.
% - beta [complex]: coefficient before t in exponential.
%
% Outputs:
% - integrated [complex]: result of Integration Par Partie.


    %% Compute integrated value
    
    % Define sum variable
    m = 1:(q+1);
    
    % Define intermediate matrix
    inm = (-1).^(m+1) .* factorial(q) ./ factorial(q+1-m) .* alpha ./ beta.^m ...
           .* ((b-k).^(q+1-m) .* exp(beta*b) - (a-k).^(q+1-m) .* exp(beta*a));
       
    % Sum over m to have result
    integrated = sum(inm);   
    
    % If beta is 0, other computation
    if beta == 0
        integrated = alpha / (q+1) * ((b-k)^(q+1) - (a-k)^(q+1));
    end
    
    
end
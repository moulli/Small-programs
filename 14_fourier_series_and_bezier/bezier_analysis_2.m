clear; close all; clc
addpath('functions')


%% Bezier from svg function

% Specify path
path_to_file = 'bezier_svg_2.svg';

% Import svg file
[B, BC, P] = bezier_from_svg(path_to_file);

% Define parameters
t = (0:0.01:1)';
T = length(B);

% Plot full drawing
figure
for p = 1:length(B)
    % Get points from function and t parameter
    Bt = B{p}(t);
    % Plot segment in 2D
    subplot(5, 1, 1:3)
    hold on
    scatter(Bt(:, 1), Bt(:, 2), 10, '.')
    % Plot x part of the segment
    subplot(5, 1, 4)
    hold on
    scatter(t+(p-1), Bt(:, 1), 10, '.')
    % Plot y part of the segment
    subplot(5, 1, 5)
    hold on
    scatter(t+(p-1), Bt(:, 2), 10, '.')
end
% 2D plot should have equal axes
subplot(5, 1, 1:3)
axis equal


%% Compute Fourier coefficients

% Number of positive coefficients
N = 25;

% Define fourier coefficients output
cn = zeros(2*N+1, 1);

% Loop over negative and positive coefficients
for n = -N:N
    % Loop over Bezier curves
    for p = 1:size(P, 1)
        % Get Bezier coefficients (for easier integration)
        bc = BC{p};
        % Loop over Bezier coefficients
        for bci = 1:size(bc, 1)
            % Integral limits ([0, 1] + Bezier curve number)
            a = p-1;
            b = p;
            % Power (varies between 0 and 3)
            q = bci-1;
            % Reposition Bezier curve so that it is integrated between 0 and 1
            k = p-1;
            % Coefficient before polynom in complex
            alpha = bc(bci, 1) + 1i*bc(bci, 2);
            % Exponential exposant for particular Fourier coefficient
            beta = -2*pi/T*n*1i;
            % Add to Fourier coefficient
            cn(n+N+1) = cn(n+N+1) + pep_ipp(a, b, q, k, alpha, beta)/T;
        end
    end        
end


%% Plot signal approximation

% Build full time vector
t_full = [];
for p = 1:length(B)
    t_full = cat(1, t_full, t+(p-1));
end

% Rebuild signal by simply multiplying exponentials with Fourier coefficients
signal = exp(2*pi/T*(-N:N)*1i.*t_full) * cn;

% Add to former plot
subplot(5, 1, 1:3)
hold on
plot(signal, 'k', 'LineWidth', 2)
axis equal
subplot(5, 1, 4)
hold on
plot(t_full, real(signal), 'k', 'LineWidth', 2)
subplot(5, 1, 5)
hold on
plot(t_full, imag(signal), 'k', 'LineWidth', 2)


%% Do a nice animation

% List of end points
end_pos = zeros(0, 2);

figure
drawnow
pause(10)
% Redefine time vector
t_full_animate = (0:0.05:T)';
% Loop over times
for t_f = 1:length(t_full_animate)
    % Define time
    ti = t_full_animate(t_f);
    % Loop over coefficients
    for n = 0:N
        % Find center, radius and angle
        if n == 0
            % Define center of figure
            C0 = [0, 0];
            % Get complex number
            cni = exp(2*pi/T*n*1i.*ti) * cn(N+1+n);
            % Get radius and angle of complex number
            r = abs(cni);
            a = angle(cni);
            % Find new center
            C = C0 + r*[cos(a), sin(a)];
            % Plot
%             viscircles(C0, r, 'LineWidth', 0.1, 'Color', 'k');
            plot([C0(1), C(1)], [C0(2), C(2)], 'k');
            hold on
        else
            % Else there is a pair of complexe exponentials
            for np = [-n, n]
                % Old center becomes new center
                C0 = C;
                % Get complex number
                cni = exp(2*pi/T*np*1i.*ti) * cn(N+1+np);
                % Get radius and angle of complex number
                r = abs(cni);
                a = angle(cni);
                % Find new center
                C = C0 + r*[cos(a), sin(a)];
                % Plot
                plot([C0(1), C(1)], [C0(2), C(2)], 'k');
                if np == 2 || np == -6
                    viscircles(C0, r, 'LineWidth', 0.75, 'Color', 'b');
                    plot([C0(1), C(1)], [C0(2), C(2)], 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
                end
            end        
        end
    end
    % Plot end point
    end_pos = cat(1, end_pos, C);
    scatter(end_pos(:, 1), end_pos(:, 2), 200, [0.635, 0.078, 0.184], '.')
    axis equal
    hold off
    drawnow
end


%% Plot animation output

figure
for p = 1:length(B)
    % Get points from function and t parameter
    Bt = B{p}(t);
    % Plot x part of the segment
    subplot(2, 1, 1)
    hold on
    scatter(t+(p-1), Bt(:, 1), 10, '.')
    % Plot y part of the segment
    subplot(2, 1, 2)
    hold on
    scatter(t+(p-1), Bt(:, 2), 10, '.')
end
% Plot new position
subplot(2, 1, 1)
plot(t_full, real(signal), 'k', 'LineWidth', 2)
scatter(t_full_animate, end_pos(:, 1), 200, [0.635, 0.078, 0.184], '.')
subplot(2, 1, 2)
plot(t_full, imag(signal), 'k', 'LineWidth', 2)
scatter(t_full_animate, end_pos(:, 2), 200, [0.635, 0.078, 0.184], '.')






    

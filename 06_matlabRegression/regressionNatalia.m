clear; close all; clc
addpath(genpath('/home/ljp/Programs/'));



%% Taking focus:

root = '/home/ljp/Science/GeoffreysComputer/Projects/RLS_Natalia/';
study = '';
date = '2018-12-12';
run = 'Run 02';

F = NT.Focus(root, study, date, run);

normalize_data = 0;



%% We can get the info of the file:

natalia = [F.dir('HDF5'), '/', erase(F.name, ' '), '.h5'];
stim = h5read(natalia, '/Data/Stimulus');
refco = h5read(natalia, '/Data/RefCoordinates');
dff = h5read(natalia, '/Data/Values');
if normalize_data == 1
    dff = (dff - mean(dff, 2)) ./ std(dff, [], 2);
end



%% Based on this info, we can obtain the data inside:

Params.period = 0;
[out, variables, perinf] = h5stimreg(natalia, Params);

% Plotting signals:
figure
hold on
for i = 1:5
    plot(variables(:, i))
end
legend('Constant signal equal to 1', 'Positive stimulus', 'Negative stimulus', 'Positive derivative of stimulus', 'Negative derivative of stimulus')
title('The five signals we regress against, all depending on initial stimulus', 'Interpreter', 'latex')
% Plotting approximation:
figure
subplot(3, 1, 1)
hold on
[~, n1] = max(out.R2score);
plot(perinf{3}(n1, :))
plot(sum(out.coef(n1, :) .* variables(:, 2:end), 2) + out.intercept(n1))
plot(perinf{2})
title('Highest R2 score signal, with regression and stimulus', 'Interpreter', 'latex')
subplot(3, 1, 2)
hold on
n2 = randperm(size(perinf{3}, 1), 1);
plot(perinf{3}(n2, :))
plot(sum(out.coef(n2, :) .* variables(:, 2:end), 2) + out.intercept(n2))
plot(perinf{2})
title('Random signal, with regression and stimulus', 'Interpreter', 'latex')
subplot(3, 1, 3)
hold on
[~, n3] = max(out.coef(:, 3));
plot(perinf{3}(n3, :))
plot(sum(out.coef(n3, :) .* variables(:, 2:end), 2) + out.intercept(n3))
plot(perinf{2})
title('Highest correlation to positive derivative of stimulus, with regression and stimulus', 'Interpreter', 'latex')



%% Plotting neurons with high F-stat and high coefficients:

%Parameters:
nneu = size(dff, 1);
q = 10;
% Coefficients:
coef = sqrt(sum((out.coef.^2), 2));
[maxcoef, indcoef] = sort(coef, 'descend');
indcK = sort(indcoef(1:ceil(nneu/q)));
% F-stat:
fstat = out.Fstat;
[maxfstat, indfstat] = sort(fstat, 'descend');
indfK = sort(indfstat(1:ceil(nneu/q)));
% Taking indices that are in both vectors:
compa = sum(indcK == indfK', 2);
indfin = indcK(compa == 1);
% Computing F-stat and coefficients for neurons kept:
indfinf = ((1 ./ fstat(indfin)) - (1 ./ max(fstat(indfin)))) ./ (1 ./ min(fstat(indfin)));
indfinc = ((1 ./ coef(indfin)) - (1 ./ max(coef(indfin)))) ./ (1 ./ min(coef(indfin)));
indind = (indfinf.*indfinc - min(indfinf.*indfinc)) ./ max(indfinf.*indfinc);


% Plot based on coefficients for each variable:
figure
for i = 1:4
    coeftemp = abs(out.coef(:, i));
    [maxcoef, indcoef] = sort(coeftemp, 'descend');
    indcK = sort(indcoef(1:ceil(nneu/q)));
    compa = sum(indcK == indfK', 2);
    indfintemp = indcK(compa == 1);
    indfinctemp = ((1 ./ coef(indfintemp)) - (1 ./ max(coef(indfintemp)))) ./ (1 ./ min(coef(indfintemp)));
    subplot(2, 2, i)
    hold on
    grid on
    axis equal
    scatter3(refco(indfintemp, 1), refco(indfintemp, 2), refco(indfintemp, 3), [], [0.9*ones(size(indfinctemp)), indfinctemp, 0.9*ones(size(indfinctemp))], '.')
    axis equal
    xlabel('x-coordinate', 'Interpreter', 'latex')
    ylabel('y-coordinate', 'Interpreter', 'latex')
    zlabel('z-coordinate', 'Interpreter', 'latex')
end
subplot(2, 2, 1)
title('Neurons with high F-stat and high absolute coefficient to positive stimulus', 'Interpreter', 'latex')
subplot(2, 2, 2)
title('Neurons with high F-stat and high absolute coefficient to negative stimulus', 'Interpreter', 'latex')
subplot(2, 2, 3)
title('Neurons with high F-stat and high absolute coefficient to positive derivative of stimulus', 'Interpreter', 'latex')
subplot(2, 2, 4)
title('Neurons with high F-stat and high absolute coefficient to negative derivative of stimulus', 'Interpreter', 'latex')


figure
subplot(1, 2, 1)
hold on
grid on
axis equal
scatter3(refco(indfin, 1), refco(indfin, 2), refco(indfin, 3), [], [indfinf, 0.9*ones(size(indfinf)), 0.9*ones(size(indfinf))], '.')
axis equal
title('Neurons with high F-stat', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')
subplot(1, 2, 2)
hold on
grid on
axis equal
scatter3(refco(indfin, 1), refco(indfin, 2), refco(indfin, 3), [], [0.9*ones(size(indfinf)), indfinc, 0.9*ones(size(indfinf))], '.')
axis equal
title('Neurons with high regression coefficients', 'Interpreter', 'latex')
xlabel('x-coordinate', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')
zlabel('z-coordinate', 'Interpreter', 'latex')











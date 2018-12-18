clear; close all; clc



%% Going to the right address:

sinstim = 0;
switch sinstim
    case 0
        cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-05-24/Run 08/Analysis/HDF5'
    case 1
        cd '/home/ljp/Science/GeoffreysComputer/RLS/Data/2018-06-14/Run 03/Analysis/HDF5'
end




%% We can get the info of the file:

hippo = 'h5hippo.h5';
stim = h5read(hippo, '/Data/Stimulus');
dff = h5read(hippo, '/Data/Values');
dffn = (dff - mean(dff, 2)) ./ std(dff, [], 2);



%% Based on this info, we can obtain the data inside:

addpath(genpath('~/Science/Hippolyte/Small-programs'))
Params.period = 1;
Params.nperiod = 10;
[out, variables, perinf] = h5stimreg(hippo, Params);
figure
subplot(2, 1, 1)
hold on
[~, n1] = max(out.R2score);
plot(perinf{3}(n1, :))
plot(sum(out.coef(n1, :) .* variables(:, 2:end), 2) + out.intercept(n1))
plot(perinf{2})
title('Highest R2 score signal, with regression and stimulus', 'Interpreter', 'latex')
subplot(2, 1, 2)
hold on
n2 = randperm(size(perinf{3}, 1), 1);
plot(perinf{3}(n2, :))
plot(sum(out.coef(n2, :) .* variables(:, 2:end), 2) + out.intercept(n2))
plot(perinf{2})
title('Random signal, with regression and stimulus', 'Interpreter', 'latex')






clear; close all; clc



%% Metropolis-Hastings algorithm:

ls = 50;
epoch = 10000;

h = rand(1, ls);
h(9) = 10;
J = 2 * randn(ls, ls); J = (J + J') ./ 2; J = J - diag(diag(J));

tot = zeros(ls, epoch);
tot(:, 1) = 2 * round(rand(ls, 1)) - 1;
pold = exp(h * tot(:, 1) + tot(:, 1)' * J * tot(:, 1));

for i = 2:epoch
    rpt = ceil(ls * rand(1, 1));
    staten = tot(:, i-1);
    staten(rpt) = (staten(rpt) == 0) .* 1;
    pnew = exp(h * staten + staten' * J * staten); 
    randu = rand(1, 1);
    if pnew / pold >= randu
        tot(:, i) = staten;
    else
        tot(:, i) = tot(:, i-1);
    end
end

figure
subplot(3, 1, 1)
plot(mean(tot, 2), 'o:')
title('Mean activity of each neuron', 'Interpreter', 'latex')
xlabel('Neuron', 'Interpreter', 'latex')
ylabel('Mean activity', 'Interpreter', 'latex')
grid on
subplot(3, 1, 2:3)
image(tot, 'CDataMapping', 'scaled')
title('Activity in time for each neuron', 'Interpreter', 'latex')
xlabel('Time', 'Interpreter', 'latex')
ylabel('Neuron', 'Interpreter', 'latex')



%% With two populations:

ls = 25;
epoch = 10000;

h = 5 * rand(1, 2 * ls);
J = 2 * randn(2 * ls, 2 * ls) - 5; J = (J + J') ./ 2; 
J(1:ls, 1:ls) = 2 * randn(ls, ls) + 5;
J((ls+1):end, (ls+1):end) = 2 * randn(ls, ls) + 5;
J = J - diag(diag(J));

tot = zeros(2 * ls, epoch);
tot(:, 1) = round(rand(2 * ls, 1));
pold = exp(h * tot(:, 1) + tot(:, 1)' * J * tot(:, 1));

for i = 2:epoch
    rpt = ceil(2 * ls * rand(1, 1));
    staten = tot(:, i-1);
    staten(rpt) = (staten(rpt) == 0) .* 1;
    pnew = exp(h * staten + staten' * J * staten); 
    randu = rand(1, 1);
    if pnew / pold >= randu
        tot(:, i) = staten;
    else
        tot(:, i) = tot(:, i-1);
    end
end

figure
subplot(6, 2, [1, 3])
image(h, 'CDataMapping', 'scaled')
colorbar
title('Solo vector', 'Interpreter', 'latex')
xlabel('Neuron', 'Interpreter', 'latex')
subplot(6, 2, [2, 4])
image(J, 'CDataMapping', 'scaled')
colorbar
title('Pairwise matrix', 'Interpreter', 'latex')
xlabel('Neuron', 'Interpreter', 'latex')
ylabel('Neuron', 'Interpreter', 'latex')
subplot(6, 2, 5:6)
plot(mean(tot, 2), 'o:')
title('Mean activity of each neuron', 'Interpreter', 'latex')
xlabel('Neuron', 'Interpreter', 'latex')
ylabel('Mean activity', 'Interpreter', 'latex')
grid on
subplot(6, 2, 7:10)
image(tot, 'CDataMapping', 'scaled')
title('Activity in time for each neuron', 'Interpreter', 'latex')
xlabel('Time', 'Interpreter', 'latex')
ylabel('Neuron', 'Interpreter', 'latex')
subplot(6, 2, 11:12)
hold on
plot(mean(tot(1:ls, :)))
plot(mean(tot(ls+1:end, :)))
legend('Population 1', 'Population 2')
title('Mean activity of both populations', 'Interpreter', 'latex')
xlabel('Time', 'Interpreter', 'latex')

        
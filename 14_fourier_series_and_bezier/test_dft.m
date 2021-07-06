clear; close all; clc

addpath(genpath('C:\Users\Hippolyte Moulle\Documents\GitHub\Small-programs'));

t = 0:0.4:100;
signal = [10 + sqrt(t).*cos(2*pi*1*t) + 2*sin(2*pi*0.3*t+pi/4);
          -4 + cos(2*pi*0.5*t+pi/12) - log(t+1).*cos(2*pi*0.01*t) + sin(2*pi*1.1*t)];
fs = 2.5;
fcut = 0.9;

[f, A, phi, low_pass] = dft(signal, fs, fcut);

figure
subplot(2, 1, 1)
hold on
for i = 1:size(A, 1)
    plot(f, A(i, :))
end
subplot(2, 1, 2)
hold on
for i = 1:size(phi, 1)
    plot(f, phi(i, :))
end

signal_rebuilt = low_pass(t);
figure
hold on
for i = 1:size(signal, 1)
    plot(signal(i, :))
    plot(signal_rebuilt(i, :))
end
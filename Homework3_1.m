%% * Homework3_1*
%% Clear workspace
close all;
clear;
clc;
%% Programmers
% Mohammad Mahdi Elyasi - 9823007
%
% Moein Nasiri - 9823093
%% Homework1
% Here we want to use diffrential equation and Z-Transform
%
% Here we declare some basic variables
fs = 10000;
f0 = 500;
w0 = 2 * pi * f0 / fs;
R = [0.8 0.9 0.99];
%%%
% Here we want to plot H for different R and see how it
% responds.
for i = 1:3

    G = (1 - R(i)) * (1 - 2 * R(i) * cos(2 * w0) + R(i) ^ 2) ^ (0.5);
    w = linspace(0, pi, 500);
    H = freqz(G, [1 -2 * R(i) * cos(w0) R(i) ^ 2], w);
    figure('Name', 'Impulse Response power 2');
    subplot(2, 1, 1)
    plot(w, abs(H) .^ 2);
    title('Absolute Impulse response power 2');
    xlabel('Frequency');
    ylabel('Amplitude');
    grid on;
    subplot(2, 1, 2)
    plot(w, angle(H) .^ 2);
    title('Phase of Impulse response power 2');
    xlabel('Frequency');
    ylabel('Phase');
    grid on;
end


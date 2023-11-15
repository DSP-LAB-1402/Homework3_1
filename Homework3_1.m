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

%% Homework2
w_2 = linspace(0, pi, 301);

for j = 1:3
    x = [1 zeros(1, 300)];
    w1 = 0;
    w2 = 0;
    y = zeros(1, 301);

    for i = 1:301
        y(i) = 2 * R(j) * cos(w0) * w1 -R(j) ^ 2 * w2 + G * x(i);
        w2 = w1;
        w1 = y(i);
    end

    figure('Name', 'Output');
    subplot(3, 1, 1)
    stem(y);
    title('Output');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    subplot(3, 1, 2)
    stem(abs(y));
    title('Absolute of Output');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    subplot(3, 1, 3)
    stem(angle(y));
    title('Phase of Output');
    xlabel('Time');
    ylabel('Phase');
    grid on;

end

for i = 1:3
    filtered_signal = filter(G, [1 -2 * R(i) .* cos(w0) R(i) .^ 2], x);
    figure('Name', 'Filtered Signal');
    subplot(3, 1, 1)
    stem(filtered_signal);
    title('Filtered Signal');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    subplot(3, 1, 2)
    stem(abs(filtered_signal));
    title('Absolute of Filtered Signal');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    subplot(3, 1, 3)
    stem(angle(filtered_signal));
    title('Phase of Filtered Signal');
    xlabel('Time');
    ylabel('Phase');
    grid on;
end

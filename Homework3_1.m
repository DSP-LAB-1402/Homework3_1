%% *Homework3_1 *
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
% In this task we want to show how impulse input responds to our filter
% by 3 different methods
%
% Firstly, we declare essential variables
for j = 1:3
    x = [1 zeros(1, 300)];
    w1 = 0;
    w2 = 0;
    y = zeros(1, 301);
    %%%
    % Here we show output by differential equation and plot it
    for i = 1:301
        y(i) = -2 * R(j) * cos(w0) * w1 -R(j) ^ 2 * w2 + G * x(i);
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

%%%
% Here we simply filtered signal by filter we created and plot it
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

%% Homework3

x = 0.1 * randn(1, 301);

for j = 1:3
    w1 = 0;
    w2 = 0;
    w21 = 0;
    w11 = 0;
    y = zeros(1, 301);
    y1 = zeros(1, 301);
    x1 = [1 zeros(1, 300)];

    for i = 1:301
        y(i) = -2 * R(j) * cos(w0) * w1 -R(j) ^ 2 * w2 + G * x(i);
        y1(i) = -2 * R(j) * cos(w0) * w11 -R(j) ^ 2 * w21 + G * x1(i);
        w2 = w1;
        w1 = y(i);
        w21 = w11;
        w11 = y1(i);
    end

    figure('Name', 'Output');
    subplot(3, 1, 1)
    stem(y);
    hold on;
    stem(y1);
    title('Output');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    legend('Noise', 'Impulse Response');
    subplot(3, 1, 2)
    stem(abs(y));
    hold on;
    stem(abs(y1));
    title('Absolute of Output');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    legend('Noise', 'Impulse Response');
    subplot(3, 1, 3)
    stem(angle(y));
    hold on;
    stem(angle(y1));
    title('Phase of Output');
    xlabel('Time');
    ylabel('Phase');
    grid on;
    legend('Noise', 'Impulse Response');
end

%% Homework4

for j = 1:3
    w1 = 0;
    w2 = 0;
    y = zeros(1, 301);

    for i = 1:301
        y(i) = -2 * R(j) * cos(w0) * w1 -R(j) ^ 2 * w2 + G * x(i);
        w2 = w1;
        w1 = y(i);
    end

    figure('Name', 'Output');
    subplot(3, 1, 1)
    stem(y);
    hold on;
    stem(x);
    title('Output');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    legend('Filtered Noise', 'Noise');
    subplot(3, 1, 2)
    stem(abs(y));
    hold on;
    stem(abs(x));
    title('Absolute of Output');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
    legend('Filtered Noise', 'Noise');
    subplot(3, 1, 3)
    stem(angle(y));
    hold on;
    stem(angle(x));
    title('Phase of Output');
    xlabel('Time');
    ylabel('Phase');
    grid on;
    legend('Filtered Noise', 'Noise');
end

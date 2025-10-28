clear;clc;

% Setting up the time vector
t = linspace(0,36,1000);

% Functions for velocity (from Part A) and altitude (antiderivative of
% velocity function, calculated using Wolfram Alpha)
vel = @(t) 2000.*log(18./(18-0.4*t))-9.81.*t;
altitude = @(t) 2000.*t - 4.905.*(t.^2) + (2000.*t - 90000).*log(45./(45-t));

% Evlauating velocity and altitude along t
v = vel(t);
alt = altitude(t);

% Plotting
tiledlayout(1,2, TileSpacing='tight', Padding='tight')

nexttile
hold on
title("Velocity vs Time (no drag)")
xlabel("time (s)")
ylabel("Velocity (m/s)")
plot(t, v, '-b')

nexttile
hold on
title("Altitude vs Time")
subtitle("MATLAB graph is without drag, hand-drawn is with drag")
xlabel("time (s)")
ylabel("Altitude (m)")
plot(t, alt, '-r')



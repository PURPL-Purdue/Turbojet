function [T8,P8,Ve] = Nozzle(T5,P5,PS0,M,N)
c = 1156; % specific heat for hot components [J/kgK]
gamma = 1.33;
eff_nozzle = 0.95; % nozzle efficiency
T8 = T5;
P8 = P5;
Ve = sqrt(2 * c * T5 * eff_nozzle * ...
    (1- (PS0 / P5)^((gamma - 1)/gamma)));       % exit velocity assuming perfectly expanded [m/s]

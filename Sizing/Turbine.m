function [T5,P5]=Turbine(T4,P4,T3,T2,FAR,M,N)
%need to add adjustments for alternate shaft speeds
c = 1156; % specific heat for hot components [J/kgK]
c_cold = 1004; % specific heat for cold components [J/kgK]
gamma = 1.33;
eff_turbine = 0.9; % turbine efficiency
T5 = (((1 + FAR)*T4*c) - ...
    (c_cold*(T3 - T2))) / ...
    ((1 + FAR) * c);                                     % stagnation temp exiting turbine [K] 
P5 = P4 * (1 - (1- T5/T4) ...
    /eff_turbine) ^ (gamma / (gamma - 1));                         % stagnation pressure exiting turbine [kPa]

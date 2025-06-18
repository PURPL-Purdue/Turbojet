function [T2,P2]=Diffuser(T0,PS0,M)
gamma = 1.4; % specific heat ratio for cold components [ambient, diffuser, compressor]
eff_diffuser = 0.95; % diffuser efficiency
T2 = T0;                                                % stagnation temp exiting diffuser [K]
P2 = PS0 * (1 + eff_diffuser*(gamma - 1) ...
    / 2 * M^2) ^ (gamma / (gamma - 1));                   % stagnation pressure exiting diffuser [kPa]

function [T3,P3]=Compressor(T2,P2,M,N)
%need to add adjustments for alternate shaft speeds
gamma = 1.4; % specific heat ratio for cold components [ambient, diffuser, compressor]
eff_compressor = 1; % compressor efficiency
CPR = 3; % compressor pressure ratio
T3 = T2 * (1 + ((1/eff_compressor) * ...
     ((CPR ^ ((gamma - 1)/gamma)) - 1)));              % stagnation temp exiting compressor [K]
P3 = CPR * P2;                             % stagnation pressure exiting compressor [kPa]

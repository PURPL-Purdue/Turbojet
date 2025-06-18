function [T4,P4,FAR]=Combuster(T3,P3)
c = 1156; % specific heat for hot components [J/kgK]]
T4 = 977; % Burner exit temperature [K]
eff_burner = 0.9; % burner efficiency
heatrxn_fuel = 50350000; % heat of combustion of jetA [J/kg]
P4 = P3;                                                % stagnation pressure exiting burner [kPa]
FAR = ((T4/T3) - 1) / (((eff_burner*heatrxn_fuel) ...
    / (c*T3)) - (T4/T3));              % fuel to air ratio
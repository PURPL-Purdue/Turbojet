%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPL
% Turbojet Sizing Code
% Created by: Josh Hyatt
%
% Program Description
% This program generates performance values for turbojet sizing
% using CFTurbo software 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION 

mach_num = 0; % ambient mach number 
gamma_cold = 1.4; % specific heat ratio for cold components [ambient, diffuser, compressor]
gamma_hot = 1.33; % specific heat ratio for hot components [burner, turbine, nozzle]
specHeat_cold = 1004; % specific heat for cold components [J/kgK]
specHeat_hot = 1156; % specific heat for hot components [J/kgK]
heatrxn_fuel = 42800000; % heat of combustion of jetA [J/kg]

eff_diffuser = 0.95; % diffuser efficiency
eff_compressor = 0.44; % compressor efficiency
eff_burner = 0.9; % burner efficiency
eff_turbine = 0.9; % turbine efficiency
eff_nozzle = 0.95; % nozzle efficiency

compr_pressure = 4; % compressor pressure ratio
burner_temp = 1300; % Burner inlet temperature [K]

stat_temp0 = 293.15; % initial static temp [K]
stat_press0 = 100; % initial static pressure [kPa]

thrust_lbf = 50; % desired thrust [lbf]
thrust_N = thrust_lbf * 4.448; % desired thrust [N]

%% CALCULATIONS 

% STAGNATION VALUES 
stagn_temp = stat_temp0 * (1 + ((gamma_cold - 1)/2 * mach_num^2));         % initial stagnation temp [K]

% DIFFUSER OUTLET / COMPRESSOR INLET (ASSUMING ISENTROPIC)
diffuser_temp = stagn_temp;                                                % stagnation temp exiting diffuser [K]
diffuser_press = stat_press0 * (1 + eff_diffuser*(gamma_cold - 1) ...
    / 2 * mach_num^2) ^ (gamma_cold / (gamma_cold - 1));                   % stagnation pressure exiting diffuser [kPa]

% COMPRESSOR OUTLET / BURNER INLET
compr_temp = diffuser_temp * (1 + ((1/eff_compressor) * ...
    ((compr_pressure ^ ((gamma_cold - 1)/gamma_cold)) - 1)));              % stagnation temp exiting compressor [K]
compr_press = compr_pressure * diffuser_press;                             % stagnation pressure exiting compressor [kPa]

% BURNER OUTLET / TURBINE INLET
burner_press = compr_press;                                                % stagnation pressure exiting burner [kPa]
fuel_ratio = ((burner_temp/compr_temp) - 1) / (((eff_burner*heatrxn_fuel) ...
    / (specHeat_hot*compr_temp)) - (burner_temp/compr_temp));              % fuel to air ratio

% TURBINE OUTLET / NOZZLE INLET
turbine_temp = (((1 + fuel_ratio)*burner_temp*specHeat_hot) - ...
    (specHeat_cold*(compr_temp - diffuser_temp))) / ...
    ((1 + fuel_ratio) * specHeat_hot);                                     % stagnation temp exiting turbine [K] 
turbine_press = burner_press * (1 - (1- turbine_temp/burner_temp) ...
    /eff_turbine) ^ (gamma_hot / (gamma_hot - 1));                         % stagnation pressure exiting turbine [kPa]

% NOZZLE OUTLET
exit_velocity = sqrt(2 * specHeat_hot * turbine_temp * eff_nozzle * ...
    (1- (stat_press0 / turbine_press)^((gamma_hot - 1)/gamma_hot)));       % exit velocity assuming perfectly expanded [m/s]

% PERFORMANCE VALUES 
specific_thrust = (1 + fuel_ratio) * exit_velocity;                        % specific thrust assuming static operation [m/s]
st_fuel = fuel_ratio / specific_thrust;                                    % Thrust Specific Fuel Consumption [kg/Ns]
eff_thermal = ((1 + fuel_ratio) * (exit_velocity^2 / 2)) ...
    / (fuel_ratio * heatrxn_fuel);                                         % Thermal Efficiency 

% MASS FLOW 
syms mass_flow
mass_flow = vpa(solve(thrust_N == specific_thrust * mass_flow,mass_flow)); % exit mass flow [kg/s]
mass_flow_air = mass_flow / (1 + fuel_ratio);                              % mass flow of air [kg/s]
mass_flow_fuel = mass_flow - mass_flow_air;                                % mass flow of fuel [kg/s]

%% RESULTS 

fprintf("\n-----------------TURBOJET PERFORMANCE PARAMETERS-----------------\n\n");
disp("Type of Fuel Used: Jet-A");
fprintf("\nDiffuser Outlet/ Compressor Inlet Temperature: %f K\n", diffuser_temp);
fprintf("Diffuser Outlet/ Compressor Inlet Pressure: %f kPa\n\n", diffuser_press);
fprintf("Compressor Outlet/ Burner Inlet Temperature: %f K\n", compr_temp);
fprintf("Compressor Outlet/ Burner Inlet Pressure: %f kPa\n\n", compr_press);
fprintf("Burner Outlet/ Turbine Inlet Temperature: %f K\n", burner_temp);
fprintf("Burner Outlet/ Turbine Inlet Pressure: %f kPa\n\n", burner_press);
fprintf("Turbine Outlet/ Nozzle Inlet Temperature: %f K\n", turbine_temp);
fprintf("Turbine Outlet/ Nozzle Inlet Pressure: %f kPa\n\n", turbine_press);
fprintf("Fuel to Air Ratio: %f\n", fuel_ratio);
fprintf("Nozzle Exit Velocity: %f m/s\n", exit_velocity);
fprintf("Nozzle Exit Specific Thrust: %f m/s\n\n", specific_thrust);
fprintf("Nozzle Exit Mass Flow: %f kg/s\n", mass_flow);
fprintf("Mass Flow of Air: %f kg/s\n", mass_flow_air);
fprintf("Mass Flow of Fuel: %f kg/s\n\n", mass_flow_fuel);
fprintf("Thrust Specific Fuel Consumption: %f kg/Ns\n", st_fuel);
fprintf("Thermal Efficiency: %f\n", eff_thermal);
fprintf("\n");


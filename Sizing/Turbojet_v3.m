
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPL
% Turbojet Sizing Code
% Created by: Josh Hyatt
% Update by: David Salama
%
% Program Description
% This program generates performance values for turbojet sizing
% using CFTurbo software 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION 
Altitude = input("what altitude are you opperating at (ft)? ");
N = input("What Percent shaft speed do you want to run at? ");


M = 0.164; % ambient mach number 
gamma_cold = 1.4; % specific heat ratio for cold components [ambient, diffuser, compressor]

heatrxn_fuel = 50350000; % heat of combustion of jetA [J/kg]

V0 = 360.32517; % Relative Velocity at the inlet
Ae = 2.5^2*pi; % %area at the exit of the nozzle

[TS0_R,PS0_A,rho] = atmosphere(Altitude,0);
TS0 = TS0_R * 5/9;
PS0 = PS0_A /20.885;

thrust_lbf = 50; % desired thrust [lbf]
thrust_N = thrust_lbf * 4.448; % desired thrust [N]

%% CALCULATIONS 

% STAGNATION VALUES 
T0 = TS0 * (1 + ((gamma_cold - 1)/2 * M^2));         % initial stagnation temp [K]

% DIFFUSER OUTLET / COMPRESSOR INLET (ASSUMING ISENTROPIC)
[T2,P2] = Diffuser(T0,PS0,M);

% COMPRESSOR OUTLET / BURNER INLET
[T3,P3] = Compressor(T2,P2,M,N);

% BURNER OUTLET / TURBINE INLET
[T4,P4,FAR] = Combuster(T3,P3);

% TURBINE OUTLET / NOZZLE INLET
[T5,P5] = Turbine(T4,P4,T3,T2,FAR,M,N);

% NOZZLE OUTLET
[T8,P8,Ve] = Nozzle(T5,P5,PS0,M,N);

% PERFORMANCE VALUES 
specific_thrust = (1 + FAR) * Ve;                        % specific thrust assuming static operation [m/s]
st_fuel = FAR / specific_thrust;                                    % Thrust Specific Fuel Consumption [kg/Ns]
eff_thermal = ((1 + FAR) * (Ve^2 / 2)) ...
    / (FAR * heatrxn_fuel);                                         % Thermal Efficiency 

% MASS FLOW 
syms mass_flow
mass_flow = vpa(solve(thrust_N == specific_thrust * mass_flow,mass_flow)); % exit mass flow [kg/s]
mass_flow_air = mass_flow / (1 + FAR);                              % mass flow of air [kg/s]
mass_flow_fuel = mass_flow - mass_flow_air;                                % mass flow of fuel [kg/s]


F = (mass_flow*Ve) - (mass_flow_air*V0); %+(P8-PS0)*Ae;

%% RESULTS 

fprintf("\n-----------------TURBOJET PERFORMANCE PARAMETERS-----------------\n\n");
disp("Type of Fuel Used: Jet-A");
fprintf("\nDiffuser Outlet/ Compressor Inlet Temperature: %f K\n", T2);
fprintf("Diffuser Outlet/ Compressor Inlet Pressure: %f kPa\n\n", P2);
fprintf("Compressor Outlet/ Burner Inlet Temperature: %f K\n", T3);
fprintf("Compressor Outlet/ Burner Inlet Pressure: %f kPa\n\n", P3);
fprintf("Burner Outlet/ Turbine Inlet Temperature: %f K\n", T4);
fprintf("Burner Outlet/ Turbine Inlet Pressure: %f kPa\n\n", P4);
fprintf("Turbine Outlet/ Nozzle Inlet Temperature: %f K\n", T5);
fprintf("Turbine Outlet/ Nozzle Inlet Pressure: %f kPa\n\n", P5);
fprintf("Fuel to Air Ratio: %f\n", FAR);
fprintf("Nozzle Exit Velocity: %f m/s\n", Ve);
fprintf("Nozzle Exit Specific Thrust: %f m/s\n\n", specific_thrust);
fprintf("Nozzle Exit Mass Flow: %f kg/s\n", mass_flow);
fprintf("Mass Flow of Air: %f kg/s\n", mass_flow_air);
fprintf("Mass Flow of Fuel: %f kg/s\n\n", mass_flow_fuel);
fprintf("Thrust Specific Fuel Consumption: %f kg/Ns\n", st_fuel);
fprintf("Thermal Efficiency: %f\n\n", eff_thermal);
fprintf("Calculated Thrust: %f lbf \n", F);
fprintf("\n");


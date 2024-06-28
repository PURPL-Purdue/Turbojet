%PURPL Turbojet Sizing
%Josh Hyatt
%5/23/2024

M0 = 0; %mach number 
gamma_cold = 1.4; %specific heat ratio for cold components
gamma_hot = 1.33; %specific heat ratio for hot components
cpcold = 1005; %specific heat for cold components [J/kgK]
cphot = 1156; %specific heat for hot components [J/kgK]
deltahb = 42800000; %heat of combustion of jetA [J/kg]

etad = 0.95; %diffuser efficiency
etan = 0.95; %nozzle efficiency
etac = 0.45; %compressor efficiency
etat = 0.9; %turbine efficiency
etab = 0.9; %burner efficiency

pic = 3; %compressor pressure ratio
t04 = 1000; %Burner temperature [K]

t0 = 293.15; %initial static temp [K]
p0 = 100; %initial static pressure [kPa]

%Stagnation values
t00 = t0 * (1+(gamma_cold-1)/2*M0^2); %initial stagnation temp [K]
p00 = p0 * (1+(gamma_cold-1)/2*M0^2)^(gamma_cold/(gamma_cold-1)); %initial stagnation pressure [kPa]

%inlet (assuming isentropic)
t02 = t00; %stagnation temp exiting inlet [K]
p02 = p0 * (1+etad*(gamma_cold-1)/2*M0^2)^(gamma_cold/(gamma_cold-1)); %stagnation pressure exiting inlet [kPa]

%compresor
t03 = t02*(1+1/(etac)*((pic)^((gamma_cold-1)/gamma_cold)-1)); %stagnation temp exiting compressor [K]
p03 = pic * p02; %stagnation pressure exiting compressor [kPa]

%combustor
p04 = p03; %stagnation pressure exiting burner [kPa]

%fuel ratio
f = ((t04/t03)-1) / ((etab*deltahb/(cphot*t03))-(t04/t03))

%Turbine
t05 = ((1+f)*t04*cphot-cpcold*(t03-t02))/((1+f)*cphot); %stagnation temp exiting turbine [K]
p05 = p04 * (1-(1-t05/t04)/etat)^(gamma_hot/(gamma_hot-1)); %stagnation pressure exiting turbine [kPa]

%calculate exit velocity assuming perfectly expanded
u9 = sqrt(2*cphot*t05*etan*(1-(p0/p05)^((gamma_hot-1)/gamma_hot))); %exit velocity [m/s]

%calculate specific thrust assuming static operation
specific_thrust = (1+f)*u9 %[m/s]

%find mass flow needed for desired thrust
thrust_lbf = 50; %desired thrust [lbf]
thrust_N = thrust_lbf * 4.448; %desired thrust [N]

syms mass_flow
mass_flow = vpa(solve(thrust_N == specific_thrust * mass_flow,mass_flow)) %exit mass flow [kg/s]

%calculte mass flow of air and fuel
mass_flow_air = mass_flow / (1+f) %mass flow of air [kg/s]
mass_flow_fuel = mass_flow - mass_flow_air %mass flow of fuel [kg/s]

t02
p02



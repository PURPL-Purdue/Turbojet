%PURPL ROTOR CALC 60K RPM
%% We have chosen single row angular contact bearings operating 
% at 60K RPM with 40 degrees celcius
% equivlent 73 AC series for skf

%% raw variables
n = 60000; % rotational speed, RPM (60K RPM)

d = 17; % inner diameter, mm (17mm)
D = 30; % outside diamter, mm (30mm)

v = 2.2; % viscousity, mm^2/s ISO VG 2 chosen

d_m = 0.5 * (d + D); % mean bearing diameter 

%% rolling friction moment calculations

% Inlet shear heating reduction factor 
phi_ish = (1) / (1 + (1.84 * (10)^(-9) ) * ( (n * d_m )^(1.28) ) * (v)^(0.64) ); 

% kinematic replenishment/stravation reduction factor

k_rs = 6 * 10^-8; %replenishment/starvations constant, grease and oil-air lubrication
k_z = 4.4; %bearing type related constant, angular contact singular row

% phi_rs
phi_rs = 1 / ...
    ( exp(1)^( k_rs * v * n * (d + D) * ...
    sqrt( (k_z / ( 2 * (D - d) ) ) ) ) );

% rolling friction variable 

% geometric constants for rolling frictional moments 
R_1 = 3.48 * 10^(-7); 
R_2 = 3.64; 
R_3 = 1.66 * 10^(-12);

%
F_g = R_3 * ((d_m)^4)*((n)^2); %

F_a = 25; % external axial load, N
F_r = 0; % external radial load, N
 
G_rr = R_1 * ( (d_m)^(1.97) ) * (F_r + F_g + R_2 * F_a)^(0.54); % variable

M_rr = phi_ish * phi_rs * G_rr * (v * n)^(0.6);  % rolling friction moment, Nmm

%% sliding frictional Moment

% sliding friction coeffecient 

mu_bl = 0.12; % constant, using movement constant
mu_ehl = 0.1; %sliding friction  coefficient in full-film conditions
% used transmission fliuds constant 

% weighting factor for the sliding friction coefficient
phi_bl = (1) / ...
    ( exp( ( (2.6 * 10^(-8) ) * ...
    ( (n * v)^(1.4) ) * d_m ) ));

mu_sl = (phi_bl * mu_bl) + (1 - phi_bl) * mu_ehl; % sliding friction coefficient 

%sliding frictional moments
s_1 = 9.85 * 10^(-3); 
s_2 = 1.55;


G_sl = s_1 * ( d_m^(0.26) ) * ( (F_r + F_g)^(4/3) + s_2 * (F_a)^(4/3) ); % sliding friction variable


% actual equation
M_sl = G_sl * mu_sl; % sliding frictional moment, Nmm

%% Frictional moement of drag

V_m = 0.0013; %drag loss factor, based on H (height of oil bath, )/ mean diameter
% based on estimated value from skf graph

i_rw = 1; %number of ball rows

K_ball = ( (i_rw * k_z * (d + D) ) / (D - d) ) * (10^(-12)); %rolling element constant

%%%

% height of oil bath
H = 1.2*d_m; % height of oil bath should be estimated as 1.2*d_m ...
% skf says to do this when H >= 1.2 * mean diameter, using outer diatmeter 
% as estimate since bearings will be fulyl submerged

%%% other constants
x = 2 * acos( ( 0.6 * d_m - H) / (0.6 * d_m) ); %

R_s = 0.36 * ( (d_m)^2 ) * (x - sin(x) ) * F_a; % 

%%% f_t fuction created

% Ai help used here, but I was able to check by having a figure ploted, it looks
% good

% Define the piecewise function
f = @(t) (t >= 0 & t <= pi) .* sin(0.5*t) + ...
         (t > pi & t < 2*pi) .* 1;
%number of periods
num_periods = 3;

% Create a vector of t values
t = linspace(0, num_periods*2*pi, 1000);

t_mod = mod(t, 2*pi);

% Evaluate the function
f_t_vals = f(t_mod);
f_t = mean(f_t_vals);
%%% final equation for frictional drag moment

M_drag = 0.4 * V_m * K_ball * ( (d_m)^2 ) * (n^2) + ...
    ( 1.093*10^(-7) ) * (n^2) * (d_m^3) * ... 
    ( ( (n * ((d_m)^2) * f_t) / v)^(-1.379) ) * R_s; %frictional drag moment Nmm



%% final equation totaling everything

M_mm = M_rr + M_sl + M_drag; %total frictional loss calcs, Nmm
M = M_mm/1000; %total frictonal loss calcs, Nm

P = (M * 2*pi * n / 60);
%%% plot 

fprintf("Rolling Frictional Moment: %.5f Nmm\n", M_rr);
fprintf("Sliding Frictional Moment: %.5f Nmm\n", M_sl);
fprintf("Frictional Moment of Drag: %.5f Nmm\n", M_drag);
fprintf("Total Frictional Moment: %.5f Nmm\n", M_mm);
fprintf("Total Frictional Moment: %.5f Nm\n", M);
fprintf("Power Loss in Watts: %.5f \n", P);
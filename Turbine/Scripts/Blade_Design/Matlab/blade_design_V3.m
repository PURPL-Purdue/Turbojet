clear;clc;clf;

%% STATOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
S_R = 65;               % Tip Radius                i| mm           5.5
S_R_LE = 1.4;           % Leading edge radius       i| mm           0.031
S_R_TE = 0.7;           % Trailing edge radius      i| mm           0.016

S_Cx = 23;              % Axial chord               i| mm           1.102
S_Ct = 1;               %  !!! Ct is no longer user specified, has a default starting value of 1 !!!
S_zeta = 4;             % Unguided turning angle    i| degrees      6.3
S_beta_IN = 0;          % Inlet blade angle         i| degrees      35
S_ep_IN = 8;            % Inlet half wedge angle    i| degrees      8
S_beta_OUT = -66;       % Exit blade angle          i| degrees      -57
S_ep_OUT = S_zeta/2;    % Exit half wedge angle     u| degrees      3.32

S_N_B = 19;             % Number of blades          i| N/A          51 
S_blade_height = 14;    % Height of blade           i| mm           N/A
stator_exclusion_factor = 0.33;

%% ROTOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
T_R = 65;               % Tip Radius                i| mm           5.5
T_R_LE = .8;            % Leading edge radius       i| mm           0.031
T_R_TE = .6;            % Trailing edge radius      i| mm           0.016

T_Cx = 20;              % Axial chord               i| mm           1.102
T_Ct = 1;               %   !!! Ct is no longer user specified, has a default starting value of 1 !!!
T_zeta = 4;             % Unguided turning angle    i| degrees      6.3
T_beta_IN = -27;        % Inlet blade angle         i| degrees      35
T_ep_IN = 10;           % Inlet half wedge angle    i| degrees      8
T_beta_OUT = -60;       % Exit blade angle          i| degrees      -57
T_ep_OUT = T_zeta/2;    % Exit half wedge angle     u| degrees      3.32

T_N_B = 31;             % Number of blades          i| N/A          51 
T_blade_height = 14;    % Height of blade           i| mm           N/A
rotor_exclusion_factor = 0.25;

%% PLOTTING OPTIONS
plot_throat = true;     % Set "true" to display the throat lines
plot_t_max  = true;     % Set "true" to display the maximum airfoil thickness lines
plot_bez_p1 = false;     % Set "true" to display the P0 -> P1 and P1 -> P2 lines
LE_align = true;        % Set "true" to align the leading edges
triangles = true;       % Set "true" to display velocity triangles

num_stators = 4;        % Number of stators to display
num_rotors  = 5;        % Number of rotors to display

% profiles_to_plot = [1,2,3];     % Which blade profiles to display                   | 1 = hub, 2 = mid, 3 = tip
% pitch_align = 2;                % Which profile to use for pitch calculations       | 1 = hub, 2 = mid, 3 = tip

%% VELOCITY TRIANGLES
velocity_triangles.V1_mag = 402;
velocity_triangles.V2_mag = 350;
velocity_triangles.V3_mag = 375;

velocity_triangles.RPM = 80000;

velocity_triangles.U_mag = velocity_triangles.RPM * 2*pi/60 * (T_R - T_blade_height/2)/1000;

velocity_triangles.a2 = 45;
velocity_triangles.B2 = -36;
velocity_triangles.B3 = -45;

velocity_triangles.W2_mag = 300;
velocity_triangles.W3_mag = 800;

% velocity_triangles.max_mag = max(v1_mag, V2_mag, V3_mag, W2_mag, W3_mag, U_mag);

%% MAIN
stator_params = blade_parameters(S_R, S_R_LE, S_R_TE, S_Cx, S_Ct, S_zeta, S_beta_IN, S_ep_IN, S_beta_OUT, S_ep_OUT, S_N_B, S_blade_height, "Stator");
rotor_params =  blade_parameters(T_R, T_R_LE, T_R_TE, T_Cx, T_Ct, T_zeta, T_beta_IN, T_ep_IN, T_beta_OUT, T_ep_OUT, T_N_B, T_blade_height, "Rotor");

stator_blade = generate_blade(stator_params, stator_exclusion_factor);
rotor_blade  = generate_blade(rotor_params, rotor_exclusion_factor);

tiledlayout(2,2, TileSpacing='tight', Padding='tight')

nexttile
title("Full Blade")
plot_set(rotor_blade, stator_blade, plot_throat, plot_t_max, plot_bez_p1, LE_align, num_stators, num_rotors, [1,2,3], 2, triangles)

nexttile
title("Hub Profiles")
plot_set(rotor_blade, stator_blade, plot_throat, plot_t_max, plot_bez_p1, LE_align, num_stators, num_rotors, [1], 1, triangles)

nexttile
title("Mid Profiles")
plot_set(rotor_blade, stator_blade, plot_throat, plot_t_max, plot_bez_p1, LE_align, num_stators, num_rotors, [2], 2, triangles)

nexttile
title("Tip Profiles")
plot_set(rotor_blade, stator_blade, plot_throat, plot_t_max, plot_bez_p1, LE_align, num_stators, num_rotors, [3], 3, triangles)

%% FUNCTIONS
% Puts all the user input into one struct for easy parameter passing
function input = blade_parameters(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, blade_height, name)
    input.R = R;
    input.R_LE = R_LE;
    input.R_TE = R_TE;
    input.Cx = Cx;
    input.Ct = Ct;
    input.zeta = zeta;
    input.beta_IN = beta_IN;
    input.ep_IN = ep_IN;
    input.beta_OUT = beta_OUT;
    input.ep_OUT = ep_OUT;
    input.N_B = N_B;
    input.blade_height = blade_height;
    input.name = name;
end
    
% Generates a blade profile with complete info set
function full_blade = generate_blade(params, exclusion_factor)
    % RADIAL EQUILIBRIUM
    params_hub = params;
    params_tip = params;
    [params_hub.beta_OUT, params_tip.beta_OUT] = rad_eq(params.beta_OUT, params.R, params.R-params.blade_height);
    [params_hub.beta_IN, params_tip.beta_IN] =   rad_eq(params.beta_IN,  params.R, params.R-params.blade_height);
    
    params_hub.R = params.R-params.blade_height;
    params.R = params.R-params.blade_height/2;

    % GENERATING GEOMETRY
    full_blade(1) = pritchard(params_hub, exclusion_factor);
    fprintf("\n" + params.name + " hub profile done\n")
    full_blade(2) = pritchard(params, exclusion_factor);
    fprintf("\n" + params.name + " mid profile done\n")
    full_blade(3) = pritchard(params_tip, exclusion_factor);
    fprintf("\n" + params.name + " tip profile done\n")
end

% Radial Equilibrium calculations
function [beta_H, beta_T] = rad_eq(theta, r_tip, r_hub)
    % info = load('e11g05n70.mat');
    % V = info.ans.res.sol.V2;
    
    V_mag = 402;
    V_theta_M = V_mag * sind(theta);
    r_mid = (r_hub + r_tip) / 2;
    
    K = V_theta_M * r_mid;
    
    V_theta_H = K/r_hub;
    V_theta_T = K/r_tip;
    
    beta_H = atand(V_theta_H/(V_mag * cosd(theta)));
    beta_T = atand(V_theta_T/(V_mag * cosd(theta)));
end
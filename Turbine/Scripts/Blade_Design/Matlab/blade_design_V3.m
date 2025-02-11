clear;clc;clf;

%% STATOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
S_R = 65;               % Tip Radius                i| mm           5.5
S_R_LE = 1.4;           % Leading edge radius       i| mm           0.031
S_R_TE = 0.7;           % Trailing edge radius      i| mm           0.016

S_ttc = 15;             % Thickness to Chord ratio  i| %            N/A
S_Cx = 17;              % Axial chord               i| mm           1.102              
S_zeta = 0.01;          % Unguided turning angle    i| degrees      6.3
S_beta_IN = 0;          % Inlet blade angle         i| degrees      35
S_ep_IN = 8;            % Inlet half wedge angle    i| degrees      8
S_beta_OUT = -66;       % Exit blade angle          i| degrees      -57
S_ep_OUT = S_zeta/2;    % Exit half wedge angle     u| degrees      3.32

S_N_B = 11;             % Number of blades          i| N/A          51 
S_blade_height = 14;    % Height of blade           i| mm           N/A
stator_exclusion_factor = 0.0;

%% ROTOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
T_R = 65;               % Tip Radius                i| mm           5.5
T_R_LE = .8;            % Leading edge radius       i| mm           0.031
T_R_TE = .6;            % Trailing edge radius      i| mm           0.016

T_ttc = 8.93;             % Thickness to Chord ratio  i| %            N/A
T_Cx = 17;              % Axial chord               i| mm           1.102
T_zeta = 0.01;          % Unguided turning angle    i| degrees      6.3
T_beta_IN = -27;        % Inlet blade angle         i| degrees      35
T_ep_IN = 10;           % Inlet half wedge angle    i| degrees      8
T_beta_OUT = -60;       % Exit blade angle          i| degrees      -57
T_ep_OUT = T_zeta/2;    % Exit half wedge angle     u| degrees      3.32

T_N_B = 17;             % Number of blades          i| N/A          51 
T_blade_height = 14;    % Height of blade           i| mm           N/A
rotor_exclusion_factor = 0.0; 

%% PLOTTING CONTROLS
plot_throat = true;     % Set "true" to display the throat lines
plot_t_max  = true;     % Set "true" to display the maximum airfoil thickness lines
plot_bez_p1 = false;     % Set "true" to display the P0 -> P1 and P1 -> P2 lines
LE_align = true;        % Set "true" to align the leading edges
triangles = true;       % Set "true" to display velocity triangles

plot_optimized_stator = true;  % Set "true" to search directory for blades of max fitness and plot those
plot_optimized_rotor = false;  % Set "true" to search directory for blades of max fitness and plot those
evo_to_search = 1;

num_stators = 4;        % Number of stators to display
num_rotors  = 5;        % Number of rotors to display

%% MAIN
if plot_optimized_stator
    stator_best_inputs = get_best(evo_to_search, "stator");
    stator_params = blade_parameters(S_R, S_R_LE, S_R_TE, stator_best_inputs(2), stator_best_inputs(1), S_zeta, S_beta_IN, S_ep_IN, S_beta_OUT, S_ep_OUT, stator_best_inputs(3), S_blade_height, "Stator");
else
    stator_params = blade_parameters(S_R, S_R_LE, S_R_TE, S_Cx, S_ttc, S_zeta, S_beta_IN, S_ep_IN, S_beta_OUT, S_ep_OUT, S_N_B, S_blade_height, "Stator");
end

if plot_optimized_rotor
    rotor_best_inputs = get_best(evo_to_search, "rotor");
    rotor_params =  blade_parameters(T_R, T_R_LE, T_R_TE, rotor_best_inputs(2), rotor_best_inputs(1), T_zeta, T_beta_IN, T_ep_IN, T_beta_OUT, T_ep_OUT, rotor_best_inputs(3), T_blade_height, "Rotor");
else
    rotor_params =  blade_parameters(T_R, T_R_LE, T_R_TE, T_Cx, T_ttc, T_zeta, T_beta_IN, T_ep_IN, T_beta_OUT, T_ep_OUT, T_N_B, T_blade_height, "Rotor");
end

stator_blade = generate_blade(stator_params, stator_exclusion_factor);
fprintf("\nOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n\n")
rotor_blade  = generate_blade(rotor_params, rotor_exclusion_factor);

tiledlayout(1,4, TileSpacing='tight', Padding='tight')

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

fprintf("\nINFO: -------------------------------\n\n")
for i = 1:3
    fprintf(stator_blade(i).parameters.name + "\n")
    fprintf("   Zweifel: " + stator_blade(i).parameters.zweifel + "\n")
    fprintf("   T_max: " + stator_blade(i).parameters.t_max + "\n")
    fprintf("   K_max_ss: " + stator_blade(i).k_max_ss + "\n")
end
fprintf("\n")
for i = 1:3
    fprintf(rotor_blade(i).parameters.name + "\n")
    fprintf("   Zweifel: " + rotor_blade(i).parameters.zweifel + "\n")
    fprintf("   T_max: " + rotor_blade(i).parameters.t_max + "\n")
    fprintf("   K_max_ss: " + rotor_blade(i).k_max_ss + "\n")
end

%% FUNCTIONS
% Puts all the user input into one struct for easy parameter passing
function input = blade_parameters(R, R_LE, R_TE, Cx, ttc, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, blade_height, name)
    input.R = R;
    input.R_LE = R_LE;
    input.R_TE = R_TE;
    input.Cx = Cx;
    input.Ct = 0;
    input.ttc = ttc;
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
    
    params_hub.name = params_hub.name + " Hub";
    params.name = params.name + " Mid";
    params_tip.name = params_tip.name + " Tip";

    params_hub.R = params.R-params.blade_height;
    params.R = params.R-params.blade_height/2;

    % GENERATING GEOMETRY
    [full_blade(1), failcode_hub] = pritchard(params_hub, exclusion_factor);
    fprintf("\n" + params_hub.name + " profile done\n")
    [full_blade(2), failcode_mid] = pritchard(params, exclusion_factor);
    fprintf("\n" + params.name + " profile done\n")
    [full_blade(3), failcode_tip] = pritchard(params_tip, exclusion_factor);
    fprintf("\n" + params_tip.name + " profile done\n")

    full_blade(1).failcode = failcode_hub;
    full_blade(2).failcode = failcode_mid;
    full_blade(3).failcode = failcode_tip;
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

function params = get_best(evolution, type)
    gen_files = dir(fullfile(type + '_optimization', type + strcat('_ev',char(string(evolution)),'gen*')));
    gens = struct('gen', {});
    for i = 1:length(gen_files)
        gens{i} = importdata(fullfile(type + '_optimization', gen_files(i).name));
    end
    
    max_fitness = 0;
    best_gen_idx = 0;
    best_bladeset_idx = 0;
    for i = 1:length(gen_files)
        gen = cell2mat(gens(i));
        for j = 1:length(gen)
            fprintf("gen: %i bladeset: %i fitness: %.3f\n", i, j, gen(j).fitness)
            fitness = gen(j).fitness;
            if fitness >= max_fitness
                max_fitness = fitness;
                best_gen_idx = i;
                best_bladeset_idx = j;
            end
        end
    end
    best_gen = cell2mat(gens(best_gen_idx));
    params = best_gen(best_bladeset_idx).inputs;
    fprintf("Best gen: %i\n", best_gen_idx)
    fprintf("Best blade: %i", best_bladeset_idx)
end
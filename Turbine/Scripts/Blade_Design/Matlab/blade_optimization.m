clear;clc;clf;

%% STATOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
S_R = 50;               % Tip Radius                i| mm           5.5
S_R_LE = 1.4;           % Leading edge radius       i| mm           0.031
S_R_TE = 0.7;           % Trailing edge radius      i| mm           0.016

% !!!! Currently being optimized for !!!!           S_ttc = 15;             % Thickness to Chord ratio  i| %            N/A
% !!!! Currently being optimiz   ed for !!!!           S_Cx = 17;              % Axial chord               i| mm           1.102              
S_zeta = 0.01;          % Unguided turning angle    i| degrees      6.3
S_beta_IN = 0;          % Inlet blade angle         i| degrees      35
S_ep_IN = 8;            % Inlet half wedge angle    i| degrees      8
S_beta_OUT = -64.48;       % Exit blade angle          i| degrees      -57
S_ep_OUT = S_zeta/2;    % Exit half wedge angle     u| degrees      3.32

% !!!! Currently being optimized for !!!!           S_N_B = 23;             % Number of blades          i| N/A          51 
S_blade_height = 10.43;    % Height of blade           i| mm           N/A
stator_exclusion_factor = 0.0;

%% ROTOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
T_R = 50;               % Tip Radius                i| mm           5.5
T_R_LE = .8;            % Leading edge radius       i| mm           0.031
T_R_TE = .6;            % Trailing edge radius      i| mm           0.016

% !!!! Currently being optimized for !!!!           T_ttc = 15;             % Thickness to Chord ratio  i| %            N/A
% !!!! Currently being optimized for !!!!           T_Cx = 15;              % Axial chord               i| mm           1.102
T_zeta = 0.01;          % Unguided turning angle    i| degrees      6.3
T_beta_IN = 23.67;        % Inlet blade angle         i| degrees      35
T_ep_IN = 10;           % Inlet half wedge angle    i| degrees      8
T_beta_OUT = -51;       % Exit blade angle          i| degrees      -57
T_ep_OUT = T_zeta/2;    % Exit half wedge angle     u| degrees      3.32

% !!!! Currently being optimized for !!!!           T_N_B = 31;             % Number of blades          i| N/A          51 
T_blade_height = 10.43;    % Height of blade           i| mm           N/A
rotor_exclusion_factor = 0.0;

%% Optimization Stuff
% Controls
blade_type = "stator";          % Which type of blade we are optimizing

evolution_number = 3;
population_size = 200;
parent_count = 2;
mutation_potency = 0.15;

% Weighting
zweifel_weight = 0.4;
curvature_weight = 0.6;

weights = [zweifel_weight, curvature_weight];

% Randomization limits
lim_ttc  = [5, 30];
lim_Cx   = [10, 25];
lim_N_B  = [13, 35];

% Get last generation, if it exists
if ~isempty(dir(fullfile(blade_type + "_optimization", blade_type + "_ev" + evolution_number + "gen*")))
    [last_gen_num, last_gen] = get_latest_gen(evolution_number, blade_type);
else
    last_gen_num = 0;
end

% Set i[ new generation
new_gen_blade = struct('inputs', {}, ...
    'profile', {}, ...
    'avg_zweifel', {}, ...
    'max_kurv', {}, ...
    'failcode', {}, ...
    'fitness', {}, ...
    'ind_name', {});

figure(1)
for i = 1:population_size
    random_input = [ ...
        rand_interval(lim_ttc), ...
        rand_interval(lim_Cx), ...
        randi(lim_N_B)];

    if last_gen_num == 0
        new_inputs = random_input;
    else
        parents = pick_parents(parent_count, last_gen);
        new_inputs = mean(vertcat(parents.genes), 1);
        new_inputs = new_inputs + ((rand(1) > 0.5)*2 - 1)*random_input.*mutation_potency;
    end

    if blade_type == "stator"
        params = blade_parameters(S_R, S_R_LE, S_R_TE, new_inputs(2), new_inputs(1), S_zeta, S_beta_IN, S_ep_IN, S_beta_OUT, S_ep_OUT, new_inputs(3), S_blade_height, blade_type);
    else
        params = blade_parameters(T_R, T_R_LE, T_R_TE, new_inputs(2), new_inputs(1), T_zeta, T_beta_IN, T_ep_IN, T_beta_OUT, T_ep_OUT, new_inputs(3), T_blade_height, blade_type);
    end
    blade = generate_blade(params, 0);

    new_gen_blade(i).inputs = new_inputs;
    new_gen_blade(i).profile = blade;
    new_gen_blade(i).avg_zweifel = (blade(1).parameters.zweifel + blade(2).parameters.zweifel + blade(3).parameters.zweifel)/3;
    new_gen_blade(i).max_kurv = max(blade(1).k_max_ss + blade(2).k_max_ss + blade(3).k_max_ss + blade(1).k_max_ps + blade(2).k_max_ps + blade(3).k_max_ps);
    if blade(1).failcode == "success!!" && blade(2).failcode == "success!!" && blade(3).failcode == "success!!"
        new_gen_blade(i).failcode = "success!!";
        plot_single(blade, true, true, false, true, 1, [1,2,3], true, blade_type)
    else
        new_gen_blade(i).failcode = "fail D:";
    end
    new_gen_blade(i).ind_name = blade_type + "_ev" + evolution_number + "gen" + (last_gen_num+1) + "no" + i;
    fprintf("\nDone generating blade number " + i + "\n")
end

for i = 1:population_size
    blade = new_gen_blade(i);
    if blade.failcode == "success!!"
        [new_gen_blade(i).fitness, minimum] = calculate_fitness(blade.avg_zweifel, blade.max_kurv, weights);
        if ~minimum
            new_gen_blade(i).fitness = 0;
            blade.failcode = "failed, didn't reach minimums";
        end
    else
        new_gen_blade(i).fitness = 0;
    end
end

gen_name = blade_type + "_ev" + char(string(evolution_number)) + "gen" + char(num2str(last_gen_num+1,'%02.f')) + ".mat";
save(blade_type + "_optimization\"+gen_name, "new_gen_blade");

load gong.mat;
soundsc(y); %Gong

%% REGULAR FUNCTIONS
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
    fprintf(params.name + " profile done\n")
    [full_blade(3), failcode_tip] = pritchard(params_tip, exclusion_factor);
    fprintf(params_tip.name + " profile done\n")

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

%% Optimization Functions (Taken from Avidh's "Optimiation Main" script)
function num = rand_interval(interval)
    num = (interval(2)-interval(1)).*rand()+interval(1);
end

function [gen_num, gen_data] = get_latest_gen(evolution, blade_type)
    gen_files = dir(fullfile(blade_type + "_optimization", blade_type + strcat('_ev',char(string(evolution)),'gen*')));
     
    gen_nums = zeros(size(gen_files));
    for i = 1:length(gen_files)
        char_name = convertStringsToChars(gen_files(i).name);
        gen_nums(i) = str2double(char_name(end-5:end-4));
    end
    
    gen_num = max(gen_nums);
    fname = blade_type + "_optimization/" + blade_type + '_ev' + num2str(evolution) + 'gen' + char(num2str(gen_num,'%02.f')) + ".mat";
    gen_data = importdata(fname);
end

function parents = pick_parents(count, gen)
    parents = struct('gen_idx', {}, 'indiv', {}, 'genes', {});
    
    for n=1:count
        weighted_probs = [gen.fitness]./sum([gen.fitness]);
        cum_weights = weighted_probs;
        for i = 1:length(weighted_probs)
            if (weighted_probs(i) > 0)
                cum_weights(i) = weighted_probs(i) + sum(weighted_probs(1:i-1));
            end
        end
    
        rand_num = rand;
    
        picked_idx = find(rand_num < cum_weights, 1);
        picked_parent = gen(picked_idx);
    
        parents(n).gen_idx = picked_idx;
        parents(n).indiv = picked_parent;
        parents(n).genes = picked_parent.inputs;

        gen(picked_idx).fitness = 0; % Remove already selected parents from consideration
    end
end

function [fitness, minimum] = calculate_fitness(zweifel, kurv, weights)

    max_kurv = 0.5;
    kurv_frac = (max_kurv-kurv)/max_kurv;
    if kurv > max_kurv
        kurv_frac = 0;
    end

    targ_zweifel = 1;
    zweifel_frac = abs(targ_zweifel - zweifel)./targ_zweifel;
    zweifel_frac = 1 - zweifel_frac;
    
    minimum = true;
    if kurv_frac < 0.2 && zweifel_frac < 0.2
        minimum = false;
    end

    norm_weights = weights./sum(weights);
    fitness = sum([zweifel_frac, kurv_frac].*norm_weights);    
end
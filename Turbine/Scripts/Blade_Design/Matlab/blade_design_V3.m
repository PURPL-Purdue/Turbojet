clear;clc;clf;

%% STATOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
S_R = 65;               % Tip Radius                i| mm           5.5
S_R_LE = 1.4;           % Leading edge radius       i| mm           0.031
S_R_TE = 0.7;           % Trailing edge radius      i| mm           0.016

S_Cx = 23;              % Axial chord               i| mm           1.102
S_Ct = 1;               % Tangential chord          i| mm           0.591
S_zeta = 4;             % Unguided turning angle    i| degrees      6.3
S_beta_IN = 0;          % Inlet blade angle         i| degrees      35
S_ep_IN = 8;            % Inlet half wedge angle    i| degrees      8
S_beta_OUT = -66;       % Exit blade angle          i| degrees      -57
S_ep_OUT = S_zeta/2;    % Exit half wedge angle     u| degrees      3.32

S_N_B = 23;             % Number of blades          i| N/A          51 
S_blade_height = 14;    % Height of blade           i| mm           N/A

%% ROTOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
T_R = 65;               % Tip Radius                i| mm           5.5
T_R_LE = 1;             % Leading edge radius       i| mm           0.031
T_R_TE = .7;            % Trailing edge radius      i| mm           0.016

T_Cx = 20;              % Axial chord               i| mm           1.102
T_Ct = 1;               % Tangential chord          i| mm           0.591
T_zeta = 4;             % Unguided turning angle    i| degrees      6.3
T_beta_IN = -20;        % Inlet blade angle         i| degrees      35
T_ep_IN = 10;           % Inlet half wedge angle    i| degrees      8
T_beta_OUT = -45;       % Exit blade angle          i| degrees      -57
T_ep_OUT = T_zeta/2;    % Exit half wedge angle     u| degrees      3.32

T_N_B = 37;             % Number of blades          i| N/A          51 
T_blade_height = 14;    % Height of blade           i| mm           N/A

%% PLOTTING OPTIONS
plot_throat = true;     % Set "true" to display the throat lines
plot_t_max  = true;     % Set "true" to display the maximum airfoil thickness lines
plot_bez_p1 = true;     % Set "true" to display the P0 -> P1 and P1 -> P2 lines

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
rotor_params =  blade_parameters(T_R, T_R_LE, T_R_TE, T_Cx, T_Ct, T_zeta, T_beta_IN, T_ep_IN, T_beta_OUT, T_ep_OUT, T_N_B, T_blade_height);
stator_params = blade_parameters(S_R, S_R_LE, S_R_TE, S_Cx, S_Ct, S_zeta, S_beta_IN, S_ep_IN, S_beta_OUT, S_ep_OUT, S_N_B, S_blade_height);

rotor_blade  = generate_blade(rotor_params);
stator_blade = generate_blade(stator_params);

plot(rotor_blade, stator_blade, plot_throat, plot_t_max, plot_bez_p1)

%% FUNCTIONS
% Puts all the user input into one struct for easy parameter passing
function input = blade_parameters(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, blade_height)
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
end

% Generates a blade profile with complete info set
function full_blade = generate_blade(params)
    % RADIAL EQUILIBRIUM
    params_hub = params;
    params_tip = params;
    [params_hub.beta_OUT, params_tip.beta_OUT] = rad_eq(params.beta_OUT, params.R, params.R-params.blade_height);
    [params_hub.beta_IN, params_tip.beta_IN] =   rad_eq(params.beta_IN,  params.R, params.R-params.blade_height);
    
    params_hub.R = params.R-params.blade_height;
    params.R = params.R-params.blade_height/2;

    % GENERATING GEOMETRY
    full_blade(1) = pritchard(params_hub);
    full_blade(2) = pritchard(params);
    full_blade(3) = pritchard(params_tip);
end

% Modifies the output matricies as required for aesthetics, then plots
function plot(rotor_blade, stator_blade, plot_throat, plot_t_max, plot_bez_p1)
    % Aligning Leading Edges
    
    % Flipping Stator upside down and translating back down
    y_flip = (rotor_blade(2).parameters.Ct + stator_blade(2).parameters.Ct)/2
    line([-2,20],[y_flip, y_flip], [stator_blade(2).R, stator_blade(2).R])
    for i = 1:3
        stator_blade(i).y_comb  = 2*y_flip - stator_blade(i).y_comb  - stator_blade(2).parameters.Ct;
        stator_blade(i).y       = 2*y_flip - stator_blade(i).y       - stator_blade(2).parameters.Ct;
        stator_blade(i).y_thicc = 2*y_flip - stator_blade(i).y_thicc - stator_blade(2).parameters.Ct;
        stator_blade(i).y_o     = 2*y_flip - stator_blade(i).y_o     - stator_blade(2).parameters.Ct;
        stator_blade(i).ss_p1y  = 2*y_flip - stator_blade(i).ss_p1y  - stator_blade(2).parameters.Ct;
        stator_blade(i).ps_p1y  = 2*y_flip - stator_blade(i).ps_p1y  - stator_blade(2).parameters.Ct;
    end
    
    for i = 1:3
        stator_blade(i).y_comb  = stator_blade(i).y_comb  - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
        stator_blade(i).y       = stator_blade(i).y       - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
        stator_blade(i).y_thicc = stator_blade(i).y_thicc - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
        stator_blade(i).y_o     = stator_blade(i).y_o     - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
        stator_blade(i).ss_p1y  = stator_blade(i).ss_p1y  - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
        stator_blade(i).ps_p1y  = stator_blade(i).ps_p1y  - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
    end
    
    %- stator_blade(2).parameters.Ct + 2*(stator_blade(i).parameters.Ct-stator_blade(2).parameters.Ct);

    % Must export at this stage in order to capture flipped stator geometry
    export(rotor_blade, stator_blade)

    % Shifting rotor to the right
    x_offset = 1.2*stator_blade(2).parameters.Cx;
    for i = 1:3
        rotor_blade(i).x_comb  = rotor_blade(i).x_comb  + x_offset;
        rotor_blade(i).x       = rotor_blade(i).x       + x_offset;
        rotor_blade(i).x_o     = rotor_blade(i).x_o     + x_offset;
        rotor_blade(i).x_thicc = rotor_blade(i).x_thicc + x_offset;
        rotor_blade(i).ss_p1x  = rotor_blade(i).ss_p1x  + x_offset;
        rotor_blade(i).ps_p1x  = rotor_blade(i).ps_p1x  + x_offset;
    end

    % PLOTTING
    [y_low, y_diff] = scale_graph(stator_blade, rotor_blade, stator_blade(2).parameters.blade_height);
    for i = 1:3
        plot_blade_V3(stator_blade(i), plot_throat, plot_t_max, plot_bez_p1)
        plot_blade_V3(rotor_blade(i), plot_throat, plot_t_max, plot_bez_p1)
    end
    
    triangle_y = y_low + 0.25*(y_diff);
    % plot_triangles(rotor_blade, stator_blade, triangle_vectors)
end
% 
function plot_triangles(rotor_blade, stator_blade, triangle_vectors)
    vel_triangle(V)
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

% Preps the XY matrix and makes it a format that Solidworks likes
function xyz = export_prep(blade, R, Ct)
    z = ones(1,length(blade.x_comb))*R;
    y = blade.y_comb-Ct/2;
    xyz = [blade.x_comb', y', z'];
    i = 1;
    while i < length(xyz)
        if xyz(i,:) == xyz(i+1,:)
            xyz(i+1,:) = [];
        elseif norm(xyz(i+1,:)-xyz(i,:)) < 0.01
            xyz(i+1,:) = [];
        else
            i = i+1;
        end
    end
    xyz = [xyz; xyz(1,:)];
end

% Writes the prepped final XY matricies to .txt files
function export(rotor_blade, stator_blade)
    % EXPORTING
    writematrix(export_prep(stator_blade(1), stator_blade(1).parameters.R, stator_blade(1).parameters.Ct), 'stator_hub.txt', 'Delimiter', 'tab');
    writematrix(export_prep(stator_blade(2), stator_blade(2).parameters.R, stator_blade(2).parameters.Ct), 'stator_mid.txt', 'Delimiter', 'tab');
    writematrix(export_prep(stator_blade(3), stator_blade(3).parameters.R, stator_blade(3).parameters.Ct), 'stator_tip.txt', 'Delimiter', 'tab');

    writematrix(export_prep(rotor_blade(1),  rotor_blade(1).parameters.R,  rotor_blade(1).parameters.Ct),  'rotor_hub.txt',  'Delimiter', 'tab');
    writematrix(export_prep(rotor_blade(2),  rotor_blade(2).parameters.R,  rotor_blade(2).parameters.Ct),  'rotor_mid.txt',  'Delimiter', 'tab');
    writematrix(export_prep(rotor_blade(3),  rotor_blade(3).parameters.R,  rotor_blade(3).parameters.Ct),  'rotor_tip.txt',  'Delimiter', 'tab');
end

% Autoscales the plot
function [y_low, y_diff] = scale_graph(stator_blade, rotor_blade, blade_height)
    figure(1)
    hold on

    x_low = min([stator_blade(1).x_comb, stator_blade(2).x_comb, stator_blade(3).x_comb]) - 2;
    x_high = max([rotor_blade(1).x_comb, rotor_blade(2).x_comb, rotor_blade(3).x_comb]) + 2;
    y_total = [stator_blade(1).y_comb, stator_blade(2).y_comb, stator_blade(3).y_comb, rotor_blade(1).y_comb, rotor_blade(2).y_comb, rotor_blade(3).y_comb];
    y_low = min(y_total) - 2;
    y_high = max(y_total) + 2;
    y_diff = y_high - y_low;
    y_low = y_low - 0.5*y_diff;
    
    xlim([x_low,x_high]);
    ylim([y_low,y_high]);
    pbaspect([x_high-x_low,y_high-y_low,blade_height]);
end
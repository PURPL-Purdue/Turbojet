clear;clc;clf;

%% STATOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
S_R = 65;             % Tip Radius                i| mm           5.5
S_R_LE = 1;           % Leading edge radius       i| mm           0.031
S_R_TE = .7;          % Trailing edge radius      i| mm           0.016

S_Cx = 23;            % Axial chord               i| mm           1.102
S_Ct = 12;            % Tangential chord          i| mm           0.591
S_zeta = 4;           % Unguided turning angle    i| degrees      6.3
S_beta_IN = 0;        % Inlet blade angle         i| degrees      35
S_ep_IN = 8;          % Inlet half wedge angle    i| degrees      8
S_beta_OUT = -50;     % Exit blade angle          i| degrees      -57
S_ep_OUT = S_zeta/2;  % Exit half wedge angle     u| degrees      3.32

S_N_B = 37;           % Number of blades          i| N/A          51 
S_blade_height = 14;  % Height of blade           i| mm           N/A

%% ROTOR INPUTS
% (The column of numbers that starts with 5.5 contain the sample values given in the Pritchard paper)
T_R = 65;             % Tip Radius                i| mm           5.5
T_R_LE = 1;           % Leading edge radius       i| mm           0.031
T_R_TE = .7;          % Trailing edge radius      i| mm           0.016

T_Cx = 23;            % Axial chord               i| mm           1.102
T_Ct = 14;            % Tangential chord          i| mm           0.591
T_zeta = 4;           % Unguided turning angle    i| degrees      6.3
T_beta_IN = 35;       % Inlet blade angle         i| degrees      35
T_ep_IN = 8;          % Inlet half wedge angle    i| degrees      8
T_beta_OUT = -50;     % Exit blade angle          i| degrees      -57
T_ep_OUT = T_zeta/2;  % Exit half wedge angle     u| degrees      3.32

T_N_B = 37;           % Number of blades          i| N/A          51 
T_blade_height = 14;  % Height of blade           i| mm           N/A

%% MAIN
rotor_params =  blade_parameters(T_R, T_R_LE, T_R_TE, T_Cx, T_Ct, T_zeta, T_beta_IN, T_ep_IN, T_beta_OUT, T_ep_OUT, T_N_B, T_blade_height);
stator_params = blade_parameters(S_R, S_R_LE, S_R_TE, S_Cx, S_Ct, S_zeta, S_beta_IN, S_ep_IN, S_beta_OUT, S_ep_OUT, S_N_B, S_blade_height);

rotor_blade  = generate_blade(rotor_params);
stator_blade = generate_blade(stator_params);

plot(rotor_blade, stator_blade)

%% FUNCTIONS
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

function full_blade = generate_blade(params)
    % RADIAL EQUILIBRIUM
    [beta_OUT_H, beta_OUT_T] = rad_eq(params.beta_OUT, params.R, params.R-params.blade_height);
    [beta_IN_H, beta_IN_T] =   rad_eq(params.beta_IN,  params.R, params.R-params.blade_height);
    
    % GENERATING GEOMETRY
    full_blade.hub = pritchard(params.R-params.blade_height,   params.R_LE, params.R_TE, params.Cx, params.Ct, params.zeta, beta_IN_H,      params.ep_IN, beta_OUT_H,      params.ep_OUT, params.N_B);
    full_blade.mid = pritchard(params.R-params.blade_height/2, params.R_LE, params.R_TE, params.Cx, params.Ct, params.zeta, params.beta_IN, params.ep_IN, params.beta_OUT, params.ep_OUT, params.N_B);
    full_blade.tip = pritchard(params.R,                       params.R_LE, params.R_TE, params.Cx, params.Ct, params.zeta, beta_IN_T,      params.ep_IN, beta_OUT_T,      params.ep_OUT, params.N_B);
    
    full_blade.params = params;
end

function plot(rotor_blade, stator_blade)
    % Flipping Stator upside down
    y_flip = (rotor_blade.params.Ct + stator_blade.params.Ct)/2;

    stator_blade.hub.y_comb = y_flip + (y_flip - stator_blade.hub.y_comb) - stator_blade.params.Ct;
    stator_blade.mid.y_comb = y_flip + (y_flip - stator_blade.mid.y_comb) - stator_blade.params.Ct;
    stator_blade.tip.y_comb = y_flip + (y_flip - stator_blade.tip.y_comb) - stator_blade.params.Ct;

    stator_blade.hub.y = y_flip + (y_flip - stator_blade.hub.y) - stator_blade.params.Ct;
    stator_blade.mid.y = y_flip + (y_flip - stator_blade.mid.y) - stator_blade.params.Ct;
    stator_blade.tip.y = y_flip + (y_flip - stator_blade.tip.y) - stator_blade.params.Ct;

    stator_blade.hub.y_thicc = y_flip + (y_flip - stator_blade.hub.y_thicc) - stator_blade.params.Ct;
    stator_blade.mid.y_thicc = y_flip + (y_flip - stator_blade.mid.y_thicc) - stator_blade.params.Ct;
    stator_blade.tip.y_thicc = y_flip + (y_flip - stator_blade.tip.y_thicc) - stator_blade.params.Ct;

    stator_blade.hub.y_o = y_flip + (y_flip - stator_blade.hub.y_o) - stator_blade.params.Ct;
    stator_blade.mid.y_o = y_flip + (y_flip - stator_blade.mid.y_o) - stator_blade.params.Ct;
    stator_blade.tip.y_o = y_flip + (y_flip - stator_blade.tip.y_o) - stator_blade.params.Ct;

    % EXPORTING
    writematrix(export_prep(stator_blade.hub, stator_blade.params.blade_height, stator_blade.params.R, stator_blade.params.Ct, 0),'stator_hub.txt','Delimiter','tab');
    writematrix(export_prep(stator_blade.mid, stator_blade.params.blade_height, stator_blade.params.R, stator_blade.params.Ct, 1),'stator_mid.txt','Delimiter','tab');
    writematrix(export_prep(stator_blade.tip, stator_blade.params.blade_height, stator_blade.params.R, stator_blade.params.Ct, 2),'stator_tip.txt','Delimiter','tab');

    writematrix(export_prep(rotor_blade.hub, rotor_blade.params.blade_height, rotor_blade.params.R, rotor_blade.params.Ct, 0),'rotor_hub.txt','Delimiter','tab');
    writematrix(export_prep(rotor_blade.mid, rotor_blade.params.blade_height, rotor_blade.params.R, rotor_blade.params.Ct, 0),'rotor_hub.txt','Delimiter','tab');
    writematrix(export_prep(rotor_blade.tip, rotor_blade.params.blade_height, rotor_blade.params.R, rotor_blade.params.Ct, 0),'rotor_hub.txt','Delimiter','tab');

    rotor_blade.hub.x_comb = rotor_blade.hub.x_comb + 1.2*stator_blade.params.Cx;
    rotor_blade.mid.x_comb = rotor_blade.mid.x_comb + 1.2*stator_blade.params.Cx;
    rotor_blade.tip.x_comb = rotor_blade.tip.x_comb + 1.2*stator_blade.params.Cx;

    % Shifting rotor to the right
    rotor_blade.hub.x = rotor_blade.hub.x + 1.2*stator_blade.params.Cx;
    rotor_blade.mid.x = rotor_blade.mid.x + 1.2*stator_blade.params.Cx;
    rotor_blade.tip.x = rotor_blade.tip.x + 1.2*stator_blade.params.Cx;
    rotor_blade(1)
    rotor_blade.hub.x_o = rotor_blade.hub.x_o + 1.2*stator_blade.params.Cx;
    rotor_blade.mid.x_o = rotor_blade.mid.x_o + 1.2*stator_blade.params.Cx;
    rotor_blade.tip.x_o = rotor_blade.tip.x_o + 1.2*stator_blade.params.Cx;

    rotor_blade.hub.x_thicc = rotor_blade.hub.x_thicc + 1.2*stator_blade.params.Cx;
    rotor_blade.mid.x_thicc = rotor_blade.mid.x_thicc + 1.2*stator_blade.params.Cx;
    rotor_blade.tip.x_thicc = rotor_blade.tip.x_thicc + 1.2*stator_blade.params.Cx;

    scale_graph(stator_blade, rotor_blade, stator_blade.params.blade_height)

    % PLOTTING
    plot_blade_V3(stator_blade.hub);
    plot_blade_V3(stator_blade.mid);
    plot_blade_V3(stator_blade.tip);
    plot_blade_V3(rotor_blade.hub);
    plot_blade_V3(rotor_blade.mid);
    plot_blade_V3(rotor_blade.tip);
end

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

function xyz = export_prep(blade, blade_height, R, Ct, multiplier)
    r_hub = R-blade_height;
    z = ones(1,length(blade.x_comb))*multiplier*(blade_height/2)+r_hub;
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

function scale_graph(stator_blade, rotor_blade, blade_height)
    figure(1)
    hold on

    x_low = min([stator_blade.hub.x_comb, stator_blade.mid.x_comb, stator_blade.tip.x_comb]) - 2;
    x_high = max([rotor_blade.hub.x_comb, rotor_blade.mid.x_comb, rotor_blade.tip.x_comb]) + 2;
    y_total = [stator_blade.hub.y_comb, stator_blade.hub.y_comb, stator_blade.hub.y_comb, rotor_blade.hub.y_comb, rotor_blade.hub.y_comb, rotor_blade.hub.y_comb];
    y_low = min(y_total) - 2;
    y_high = max(y_total) + 2;
    
    xlim([x_low,x_high]);
    ylim([y_low,y_high]);
    pbaspect([x_high-x_low,y_high-y_low,blade_height]);
end

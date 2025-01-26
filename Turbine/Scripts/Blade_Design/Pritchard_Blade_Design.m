clear;clc;clf;

%% Inputs (the column of numbers that starts with 5.5 is the sample values given in the Pritchard paper)
R = 65;             % Tip Radius                i| mm           5.5
R_LE = 1;           % Leading edge radius       i| mm           0.031
R_TE = .7;          % Trailing edge radius      i| mm           0.016

Cx = 23;            % Axial chord               i| mm           1.102
Ct = 18;            % Tangential chord          i| mm           0.591
zeta = 4;           % Unguided turning angle    i| degrees      6.3
beta_IN = 12;       % Inlet blade angle         i| degrees      35
ep_IN = 8;          % Inlet half wedge angle    i| degrees      8
beta_OUT = -50;     % Exit blade angle          i| degrees      -57
ep_OUT = zeta/2;    % Exit half wedge angle     u| degrees      3.32

N_B = 17;           % Number of blades          i| N/A          51 

blade_height = 14;
num_displayed = 1;  % Number of blades to display on the graph

%% MAIN

% [beta_OUT_H, beta_OUT_T] = rad_eq(beta_OUT, R, R-blade_height)
beta_OUT_H = beta_OUT-5;
beta_OUT_T = beta_OUT+5;

% Generating the blade X and Y matricies
hub = pritchard(R-blade_height, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT_H, ep_OUT, N_B, num_displayed);
mid = pritchard(R-blade_height/2, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, num_displayed);
tip = pritchard(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT_T, ep_OUT, N_B, num_displayed);

% Autoscaling Graph (We're (not that) cool like that) *side note autoscaling doesn't
% actually work, if someone wants to fix it and shoot me a dm on discord
% that would be really sick :D
figure(1)
hold on
x_low = -0.25;
x_high = Cx + R_LE + R_TE + 0.25;
y_low = -0.25;
y_high = (num_displayed) * 2*pi*R / N_B + 0.25;
xlim([x_low,x_high]);
ylim([y_low,y_high]);
pbaspect([x_high-x_low,y_high-y_low,blade_height]);

% Blade Plotting
plot_blade_set(hub, 0, blade_height, R);
plot_blade_set(mid, 1, blade_height, R);
plot_blade_set(tip, 2, blade_height, R);

% plot_verticals(hub, mid, tip, blade_height, R); % Plots the vertical
% lines that makes the blade look cooler in 3D. I've commented it out for
% now because using it means it takes longer to graph

% Modifies X and Y matricies for easy exporting to Solidworks, and writes
% them to .txt files for use in the "Curve through XYZ" function in
% Solidworks
writematrix(export_prep(hub(1), blade_height, R, Ct, 0),'hub.txt','Delimiter','tab');
writematrix(export_prep(mid(1), blade_height, R, Ct, 1),'mid.txt','Delimiter','tab');
writematrix(export_prep(tip(1), blade_height, R, Ct, 2),'tip.txt','Delimiter','tab');

%% MAIN FUNCTIONS
function [beta_OUT_H, beta_OUT_T] = rad_eq(alpha, r_tip, r_hub)
    % info = load('e11g05n70.mat');
    % V = info.ans.res.sol.V2;
    
    V_mag = 402;
    V_theta_M = V_mag * sind(alpha);
    r_mid = (r_hub + r_tip) / 2;
    
    K = V_theta_M * r_mid;
    
    V_theta_H = K/r_hub;
    V_theta_T = K/r_tip;
    
    beta_OUT_H = atand(V_theta_H/(V_mag * cosd(alpha)));
    beta_OUT_T = atand(V_theta_T/(V_mag * sind(alpha)));
end

function blade_set = pritchard(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, num_displayed)
    o = 2*pi*R/N_B*cos(beta_OUT*pi/180)-2*R_TE;
    pitch = 2*pi*R / N_B;
    beta_IN = deg2rad(beta_IN);
    beta_OUT = deg2rad(beta_OUT);
    ep_IN = deg2rad(ep_IN);
    ep_OUT = deg2rad(ep_OUT);
    zeta = deg2rad(zeta);    

    pts = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o);
    for  i = 1:num_displayed
        blade = bezier_curve_generator(pts.x_coords, pts.y_coords, pts.betas, o, (i-1)*pitch);
        blade = arc_generator(blade.x, blade.y, blade.betas, R_LE, R_TE, pts.R_new, blade);
        [blade.x_comb, blade.y_comb] = combiner(blade);     % Combine separate sections into big x and y vectors (good for solidworks exporting later too)
        blades(i) = blade;
    end
    blade_set = blades;
end

function plot_blade_set(blade_set, height_multiplier, blade_height, R)
    r_hub = R - blade_height;
    half_height = blade_height/2;
    for i = 1:length(blade_set)
        plot3(blade_set(i).x_comb,blade_set(i).y_comb, ones(1, length(blade_set(1).x_comb))*height_multiplier*half_height+r_hub, '-k')               % Curves!
        scatter3(blade_set(i).x, blade_set(i).y, ones(1, 6)*height_multiplier*half_height+r_hub, 4, 'b', "filled")         % Points!!
        plot3(blade_set(i).x_o, blade_set(i).y_o, ones(1, length(blade_set(1).x_o))*height_multiplier*half_height+r_hub, '--r')                   % Throat!!!
    end
end

function plot_verticals(bladeset_hub, bladeset_mid, bladeset_tip, blade_height, R)
    loft(bladeset_hub, bladeset_mid, blade_height, R, 1)
    loft(bladeset_mid, bladeset_tip, blade_height, R, 2)
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

%% HELPER FUNCTIONS
function loft(bladeset_lower, bladeset_upper, blade_height, R, multiplier)
    r_hub = R - blade_height;
    vert_xy = 1:10:length(bladeset_lower(1).x_comb);
    for i = 1:length(bladeset_lower)
        for j = vert_xy
            x = linspace(bladeset_lower(i).x_comb(j), bladeset_upper(i).x_comb(j), 10);
            y= linspace(bladeset_lower(i).y_comb(j), bladeset_upper(i).y_comb(j), 10);
            z = linspace((multiplier-1)*blade_height/2+r_hub, multiplier*blade_height/2+r_hub, 10);
            plot3(x, y, z, '-b')
        end
    end
end

function pts = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o)
    cont = 1;
    while cont
        %% Points 1-5 -> X,Y Coordinates and Tangent Angles
        % POINT ONE - Suction Surface Trailing Edge Tangency Point
        beta_1 = beta_OUT - ep_OUT;
        x1 = Cx - R_TE * (1 + sin(beta_1));
        y1 = R_TE * cos(beta_1);
        
        % POINT TWO - Suction Surface Throat Point
        beta_2 = beta_OUT - ep_OUT + zeta;
        x2 = Cx - R_TE + (o + R_TE) * sin(beta_2);
        y2 = 2*pi*R/N_B - (o + R_TE) * cos(beta_2);
        
        % POINT THREE - Suction Surface Leading Edge Tangency Point
        beta_3 = beta_IN + ep_IN;
        x3 = R_LE * (1 - sin(beta_3));
        y3 = Ct + R_LE * cos(beta_3);
        
        % POINT FOUR - Pressure Surface Leading Edge Tangency Point
        beta_4 = beta_IN - ep_IN;
        x4 = R_LE * (1 + sin(beta_4));
        y4 = Ct - R_LE * cos(beta_4);
        
        % POINT FIVE - Pressure Surface Trailing Edge Tangency Point
        beta_5 = beta_OUT + ep_OUT;
        x5 = Cx - R_TE * (1 - sin(beta_5));
        y5 = -R_TE * cos(beta_5);

        %% Iteration on ep_OUT, following method by Mansour
        syms x_perp;
        % Symbolic equation for the orthographic line at Point 1
        p1_slope = -1/tan(beta_1);
        x1_perpendicular_sym = p1_slope*(x_perp-x1)+y1;

        % Symbolic equation for the orthographic line at Point 2
        o_slope = -1/tan(beta_2);
        x2_perpendicular_sym = o_slope*(x_perp-x2)+y2;
    
        % Intersection of the two orthographic lines
        x01 = double(vpasolve(x1_perpendicular_sym == x2_perpendicular_sym, x_perp));
        y01 = double(subs(x1_perpendicular_sym, x_perp, x01));

        % Iterates on ep_OUT
        R_01 = sqrt((x1-x01)^2 + (y1-y01)^2);
        yy2 = y01 + sqrt(R_01^2 - (x2-x01)^2);
        delta = abs(yy2 - y2);
        if delta > 0.0000001
            ep_OUT = ep_OUT * (y2/yy2);
        else
            cont = 0;
        end
    end
    %% Struct Outputs
    pts.betas = [beta_1, beta_2, beta_3, beta_4, beta_5];
    pts.x_coords = [x1, x2, x3, x4, x5, x01];
    pts.y_coords = [y1, y2, y3, y4, y5, y01];
    pts.R_new = R_01;
end

function blade = bezier_curve_generator(x, y, betas, o, offset)
    y = y + offset;
    %% BEZIERS!?
    t = 0:0.001:1;
    b0 = (1-t).^2;
    b1 = 2*t.*(1-t);
    b2 = t.^2;

    % Tangent Line Equations at Points 2-5
    syms x_perp
    % Slopes
    p2_slope = tan(betas(2));
    p3_slope = tan(betas(3));
    p4_slope = tan(betas(4));
    p5_slope = tan(betas(5));
    % Tangent Line Equations
    x2_perpendicular_sym = p2_slope*(x_perp-x(2))+y(2);
    x3_perpendicular_sym = p3_slope*(x_perp-x(3))+y(3);
    x4_perpendicular_sym = p4_slope*(x_perp-x(4))+y(4);
    x5_perpendicular_sym = p5_slope*(x_perp-x(5))+y(5);

    % Intersection of tangent lines to find P1's
    ss_p1x = double(vpasolve(x2_perpendicular_sym == x3_perpendicular_sym, x_perp));
    ss_p1y = double(subs(x2_perpendicular_sym, x_perp, ss_p1x));
    ps_p1x = double(vpasolve(x4_perpendicular_sym == x5_perpendicular_sym, x_perp));
    ps_p1y = double(subs(x4_perpendicular_sym, x_perp, ps_p1x));

    % Generating X and Y Coordinates for the Bezier
    ss_x_bez = x(2)*b0 + ss_p1x*b1 + x(3)*b2;
    ss_y_bez = y(2)*b0 + ss_p1y*b1 + y(3)*b2;
    ps_x_bez = x(4)*b0 + ps_p1x*b1 + x(5)*b2;
    ps_y_bez = y(4)*b0 + ps_p1y*b1 + y(5)*b2;

    %% Throat Plotting
    o_slope = -1/tan(betas(2));                                 % Gets perpendicular slope
    x_o = linspace(x(2), x(2) + o*cos(atan(o_slope)),1000);     % X-coords for throat plot
    y_o = o_slope.*(x_o-x(2))+y(2);                             % Y-coords for throat plot

    %% Adding Points 1-5 and betas to Blade Struct
    blade.x_suction = ss_x_bez;
    blade.y_suction = ss_y_bez;
    blade.x_pressure = ps_x_bez;
    blade.y_pressure = ps_y_bez;
    blade.x_o = x_o;
    blade.y_o = y_o;
    blade.x = x;
    blade.y = y;
    blade.betas = betas;
end

function new_blade = arc_generator(x, y, betas, R_LE, R_TE, R_new, blade)
    new_blade = blade;
    [new_blade.LEx, new_blade.LEy] = arc(x(3) + R_LE * cos(pi/2 - betas(3)), y(3) - R_LE * sin(pi/2 - betas(3)), R_LE, x(3), x(4), y(3), y(4), 1);
    [new_blade.TEx, new_blade.TEy] = arc(x(1) + R_TE * cos(pi/2 - betas(1)), y(1) - R_TE * sin(pi/2 - betas(1)), R_TE, x(5), x(1), y(5), y(1), 1);
    [new_blade.UGx, new_blade.UGy] = arc(x(6), y(6), R_new, x(1), x(2), y(1), y(2), 20);
end

function [x_arc, y_arc] = arc(x,y,r, x1, x2, y1, y2, multiplier)
    start = atan2(y1-y,x1-x);
    stop = atan2(y2-y,x2-x);
    if start > stop
        stop = stop + 2 * pi;
    end

    theta = start:pi/(1000*multiplier):stop;
    x_arc = r * cos(theta) + x;
    y_arc = r * sin(theta) + y;
end

function [x_comb, y_comb] = combiner(blade)
    x_comb = [blade.x_suction, blade.LEx, blade.x_pressure, blade.TEx, blade.UGx];
    y_comb = [blade.y_suction, blade.LEy, blade.y_pressure, blade.TEy, blade.UGy];
end

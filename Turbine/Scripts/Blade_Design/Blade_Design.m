clear;clc;clf;

%% Numbers currently set are the example numbers given in the Pritchard paper, and the output matches their example output.
%% INPUTS
R = 5.5;            % Airfoil radius            i| mm           5.5
R_LE = 0.031;       % Leading edge radius       i| mm           0.031
R_TE = 0.016;       % Trailing edge radius      i| mm           0.016

Cx = 1.102;         % Axial chord               i| mm           1.102
Ct = 0.591;         % Tangential chord          i| mm           0.591
zeta = 6.3;         % Unguided turning angle    i| degrees      6.3
beta_IN = 35;       % Inlet blade angle         i| degrees      35
ep_IN = 8;          % Inlet half wedge angle    i| degrees      8
beta_OUT = -57;     % Exit blade angle          i| degrees      -57
ep_OUT = 3.32;      % Exit half wedge angle     u| degrees      3.32

N_B = 51;           % Number of blades          i| N/A          51
o = 0.337;          % Throat distance           i| mm           0.337

x_offset = 0;       % X-offset for Blade 2      i| mm
y_offset = 0.5;     % Y-offset for Blade 2      i| mm

beta_IN = deg2rad(beta_IN);
beta_OUT = deg2rad(beta_OUT);
ep_IN = deg2rad(ep_IN);
ep_OUT = deg2rad(ep_OUT);
zeta = deg2rad(zeta);

%% DOING STUFF FR
[x_coords, y_coords, betas] = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o);

blade_1 = polynomial_shit_also_blade_generator_eyes_emoji(x_coords, y_coords, betas, o, 0);
blade_2 = polynomial_shit_also_blade_generator_eyes_emoji(x_coords, y_coords, betas, o, blade_1.offset);

plot_blade(blade_1, betas, R, R_LE, R_TE);
plot_blade(blade_2, betas, R, R_LE, R_TE);

function [x_coords, y_coords, betas] = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o)
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

    % OUTPUT MATRICIES
    betas = [beta_1, beta_2, beta_3, beta_4, beta_5];
    x_coords = [x1, x2, x3, x4, x5];
    y_coords = [y1, y2, y3, y4, y5];
end
function blade = polynomial_shit_also_blade_generator_eyes_emoji(x, y, betas, o, offset)
    y = y + offset;
    %% SUCTION SIDE POLYNOMIAL
    % Coefficients (From Prichard)
    d = (tan(betas(3))+tan(betas(2)))/(x(3)-x(2))^2 - 2*(y(3)-y(2))/(x(3)-x(2))^3;
    c = (y(3)-y(2))/(x(3)-x(2))^2 - tan(betas(2))/(x(3)-x(2)) - d*(x(3) + 2*x(2));
    b = tan(betas(2)) - 2*c*x(2) - 3*d*x(2)^2;
    a = y(2) - b*x(2) - c*x(2)^2 - d*x(2)^3;
    % Polynomial
    x_suction = linspace(x(3),x(2),100);
    y_suction = a + b.*x_suction + c*x_suction.^2 + d*x_suction.^3;

    blade.x_suction = x_suction;
    blade.y_suction = y_suction;

    %% Throat Plotting Stuff
    % Symbolic Polynomial Stuff
    syms xs;
    poly_suction = a + b*xs + c*xs^2 + d*xs^3;
    poly_suction_1deriv = diff(poly_suction);               % Finds derivative of suction side polynomial
    p2_deriv = subs(poly_suction_1deriv,xs,x(2));           % Evaluates derivative at Point 2
    o_slope = -1/p2_deriv;                                  % Gets perpendicular slope
    x_o = linspace(x(2), x(2) + o*cos(atan(o_slope)),50);   % X-coords for throat plot
    y_o = o_slope.*(x_o-x(2))+y(2);                         % Y-coords for throat plot

    blade.x_o = x_o;
    blade.y_o = y_o;

    %% PRESSURE SIDE POLYNOMIAL
    % Coefficients (From Prichard)
    d = (tan(betas(4))+tan(betas(5)))/(x(4)-x(5))^2 - 2*(y(4)-y(5))/(x(4)-x(5))^3;
    c = (y(4)-y(5))/(x(4)-x(5))^2 - tan(betas(5))/(x(4)-x(5)) - d*(x(4) + 2*x(5));
    b = tan(betas(5)) - 2*c*x(5) - 3*d*x(5)^2;
    a = y(5) - b*x(5) - c*x(5)^2 - d*x(5)^3;
    % Polynomial
    x_pressure = linspace(x(4),x(5),100);
    y_pressure = a + b.*x_pressure + c*x_pressure.^2 + d*x_pressure.^3;

    blade.x_pressure = x_pressure;
    blade.y_pressure = y_pressure;

    %% Get Y Offset for Second Blade
    % Deep symbolic analysis of primary texts with Rebecca Shaw
    syms xp;
    sym_poly_p = a + b*xp + c*xp^2 + d*xp^3;
    lower_y = subs(sym_poly_p, xp, x(2) + o*cos(atan(o_slope)) );

    blade.offset = y_o(end)-lower_y;
    
    %% Adding Points 1-5 to Blade Struct
    blade.points_x = x;
    blade.points_y = y;
end
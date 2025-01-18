clear;clc;clf;

%% Numbers set are the example numbers given in the Pritchard paper, and the output matches their example output.
%% Inputs
R = 5.5;            % Disk radius               i| mm           5.5
R_LE = 0.031;       % Leading edge radius       i| mm           0.031
R_TE = 0.016;       % Trailing edge radius      i| mm           0.016

Cx = 1.102;         % Axial chord               i| mm           1.102
Ct = 0.591;         % Tangential chord          i| mm           0.591
zeta = 6.3;         % Unguided turning angle    i| degrees      6.3
beta_IN = 35;       % Inlet blade angle         i| degrees      35
ep_IN = 8;          % Inlet half wedge angle    i| degrees      8
beta_OUT = -57;     % Exit blade angle          i| degrees      -57
ep_OUT = zeta/2;    % Exit half wedge angle     u| degrees      3.32

N_B = 51;           % Number of blades          i| N/A          51
o = 0.337;          % Throat distance           i| mm           0.337

num_displayed = 2;  % Number of blades to display on the graph

beta_IN = deg2rad(beta_IN);
beta_OUT = deg2rad(beta_OUT);
ep_IN = deg2rad(ep_IN);
ep_OUT = deg2rad(ep_OUT);
zeta = deg2rad(zeta);
pitch = 2*pi*R / N_B;

%% Calculating Points and Generating Polynomials
pts = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o);

for  i = 1:num_displayed
    blades(i) = polynomial_shit_also_blade_generator_eyes_emoji(pts.x_coords, pts.y_coords, pts.betas, o, (i-1)*pitch);
end

%% Plot Setup
figure(1)
hold on
% Autoscaling Graph (We're cool like that)
x_low = -0.25;
x_high = Cx + R_LE + R_TE + 0.25;
y_low = -0.25;
y_high = (num_displayed-1) * pitch + blades(1).vert + 0.25;
xlim([x_low,x_high]);
ylim([y_low,y_high]);
pbaspect([x_high-x_low,y_high-y_low,1]);

%% Plot Blades
for i = 1:num_displayed
    plot_blade(blades(i), pts.betas, R, R_LE, R_TE, pts.R_new);
end

%% FUNCTIONS
function pts = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o)
    cont = 1;
    while cont
        %% Finding X,Y Coordinates and Tangent Angles
        % POINT ONE - Suction Surface Trailing Edge Tangency Point
        beta_1 = beta_OUT - ep_OUT;
        x1 = Cx - R_TE * (1 + sin(beta_1));
        y1 = R_TE * cos(beta_1);
        
        % POINT TWO - Suction Surface Throat Point
        beta_2 = beta_OUT - ep_OUT + zeta;
        x2 = Cx - R_TE + (o + R_TE) * sin(beta_2);
        y2 = 2*pi*R/N_B - (o + R_TE) * cos(beta_2);
        
        % POINT THREE - Suction Surface Leading Edge Tangency Point
        beta_3 = beta_IN + ep_IN; %******
        x3 = R_LE * (1 - sin(beta_3));
        y3 = Ct + R_LE * cos(beta_3);
        
        % POINT FOUR - Pressure Surface Leading Edge Tangency Point
        beta_4 = beta_IN - ep_IN;
        x4 = R_LE * (1 + sin(beta_4));
        y4 = Ct - R_LE * cos(beta_4);
        
        % POINT FIVE - Pressure Surface Trailing Edge Tangency Point
        beta_5 = beta_OUT + ep_OUT; %******
        x5 = Cx - R_TE * (1 - sin(beta_5));
        y5 = -R_TE * cos(beta_5);

        %% Finding orthographic lines from Points 1 and 2, as described in the Mansour paper (with some modifications)
        % SUCTION SIDE POLYNOMIAL
        % Coefficients (From Prichard)
        d = (tan(beta_3)+tan(beta_2))/(x3-x2)^2 - 2*(y3-y2)/(x3-x2)^3;
        c = (y3-y2)/(x3-x2)^2 - tan(beta_2)/(x3-x2) - d*(x3 + 2*x2);
        b = tan(beta_2) - 2*c*x2 - 3*d*x2^2;
        a = y2 - b*x2 - c*x2^2 - d*x2^3;
    
        syms xs;
        poly_suction = a + b*xs + c*xs^2 + d*xs^3;
        poly_suction_1deriv = diff(poly_suction);   % Symbolic cubic polynomial for the suction side
    
        syms x_perp;
        p2_deriv = subs(poly_suction_1deriv,xs,x2);
        o_slope = double(-1/p2_deriv);
        x2_perpendicular_sym = o_slope*(x_perp-x2)+y2;  % Symbolic equation for the orthographic line at Point 2
    
        p1_slope = -1/tan(beta_1);
        x1_perpendicular_sym = p1_slope*(x_perp-x1)+y1; % Symbolic equation for the orthographic line at Point 2

        % Intersection of the two orthographic lines
        x01 = double(vpasolve(x1_perpendicular_sym == x2_perpendicular_sym, x_perp));
        y01 = double(subs(x1_perpendicular_sym, x_perp, x01));

        % Iterates on ep_OUT
        R_01 = sqrt((x1-x01)^2 + (y1-y01)^2);
        yy2 = y01 + sqrt(R_01^2 - (x2-x01)^2);
        xx2 = x01 + sqrt(R_01^2 - (y2-y01)^2);
        delta = double(abs(yy2 - y2));
        if delta > 0.0000001
            ep_OUT = ep_OUT * (y2/yy2);
        else
            cont = 0;
        end
    end

    % OUTPUT STRUCT
    pts.betas = [beta_1, beta_2, beta_3, beta_4, beta_5];
    pts.x_coords = [x1, x2, x3, x4, x5, x01];
    pts.y_coords = [y1, y2, y3, y4, y5, y01];
    pts.R_new = R_01;
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
    blade.vert = max(y_suction) - y(5);

    %% Throat Plotting
    % Symbolic Polynomial Stuff
    syms xs;
    poly_suction = a + b*xs + c*xs^2 + d*xs^3;
    poly_suction_1deriv = diff(poly_suction);               % Finds deriva tive of suction side polynomial
    p2_deriv = double(subs(poly_suction_1deriv,xs,x(2)));   % Evaluates derivative at Point 2
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
    
    %% Adding Points 1-5 to Blade Struct
    blade.points_x = x;
    blade.points_y = y;
end

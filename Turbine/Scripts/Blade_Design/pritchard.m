%% GENERATES BLADE GEOMETRY
function blade = pritchard(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B)
    beta_IN = deg2rad(beta_IN);
    beta_OUT = deg2rad(beta_OUT);
    ep_IN = deg2rad(ep_IN);
    ep_OUT = deg2rad(ep_OUT);
    zeta = deg2rad(zeta);   

    o = 2*pi*R/N_B*cos(beta_OUT)-2*R_TE;

    pts = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o);
    blade = generate_curves(pts.x_coords, pts.y_coords, pts.betas, R_LE, R_TE, pts.R_new, o);
    [blade.x_comb, blade.y_comb, blade.x_ss_comb, blade.y_ss_comb] = combiner(blade);     % Combine separate sections into big x and y vectors (good for solidworks exporting later too)

    [blade.t_max, pos_ss, pos_ps] = max_t(blade);
    blade.x_thicc = linspace(blade.x_pressure(pos_ps), blade.x_ss_comb(pos_ss), 20);
    blade.y_thicc = linspace(blade.y_pressure(pos_ps), blade.y_ss_comb(pos_ss), 20);
    blade.R = R;
    fprintf(">:D\n")
end

%% HELPER FUNCTIONS
% Generates bezier curves for suction and pressure side
function [x_suction, y_suction, x_pressure, y_pressure] = generate_bezier(x, y, betas, UGx_length)
    %% BEZIERS!?

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

    % Bezier for suction side
    t = linspace(0,1,1000);
    b0 = (1-t).^2;
    b1 = 2*t.*(1-t);
    b2 = t.^2;
    x_suction = x(2)*b0 + ss_p1x*b1 + x(3)*b2;
    y_suction = y(2)*b0 + ss_p1y*b1 + y(3)*b2;
    % Bezier for pressure side (number of points matches number of combined
    % points for pressure side bezier + UGT arc
    t = linspace(0,1,1000 + UGx_length);
    b0 = (1-t).^2;
    b1 = 2*t.*(1-t);
    b2 = t.^2;
    x_pressure = x(4)*b0 + ps_p1x*b1 + x(5)*b2;
    y_pressure = y(4)*b0 + ps_p1y*b1 + y(5)*b2;
end

% Generates the throat line
function [x_o, y_o] = generate_throat(x, y, betas, o)
    %% Throat Plotting
    o_slope = -1/tan(betas(2));                                 % Gets perpendicular slope
    x_o = linspace(x(2), x(2) + o*cos(atan(o_slope)),1000);     % X-coords for throat plot
    y_o = o_slope.*(x_o-x(2))+y(2);                             % Y-coords for throat plot
end

% Generates the LE, TE, UGT circular arcs
function blade = generate_arcs(x, y, betas, R_LE, R_TE, R_new)
    [blade.LEx, blade.LEy] = arc(x(3) + R_LE * cos(pi/2 - betas(3)), y(3) - R_LE * sin(pi/2 - betas(3)), R_LE, x(3), x(4), y(3), y(4), 1);
    [blade.TEx, blade.TEy] = arc(x(1) + R_TE * cos(pi/2 - betas(1)), y(1) - R_TE * sin(pi/2 - betas(1)), R_TE, x(5), x(1), y(5), y(1), 1);
    [blade.UGx, blade.UGy] = arc(x(6), y(6), R_new, x(1), x(2), y(1), y(2), 20);
    blade.x = x;
    blade.y = y;
    blade.betas = betas;
end

% Combines all the generation functions
function blade = generate_curves(x, y, betas, R_LE, R_TE, R_new, o)
    blade = generate_arcs(x, y, betas, R_LE, R_TE, R_new);
    [blade.x_suction, blade.y_suction, blade.x_pressure, blade.y_pressure] = generate_bezier(x, y, betas, length(blade.UGx));
    [blade.x_o, blade.y_o] = generate_throat(x, y, betas, o);
end

% Generates Points 1-5 and their betas
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

% Calculates max blade thickness and its location
function [t_max, pos_ss, pos_ps] = max_t(blade)
    thiccnesses = zeros(1, length(blade.x_pressure));
    for i = 1:length(blade.x_pressure)
        thiccnesses(i) = min(vecnorm([blade.x_ss_comb(i); blade.y_ss_comb(i)] - [blade.x_pressure; blade.y_pressure]));
    end
    [t_max, pos_ss] = max(thiccnesses);
    [~, pos_ps] = min(vecnorm([blade.x_ss_comb(pos_ss); blade.y_ss_comb(pos_ss)] - [blade.x_pressure; blade.y_pressure]));
    % fprintf("t_max: %f\n", t_max)
end

% Helper function for generating arc coordinates
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

% Combines XY matricies into one big XY matrix
function [x_comb, y_comb, x_ss_comb, y_ss_comb] = combiner(blade)
    x_comb = [blade.x_suction, blade.LEx, blade.x_pressure, blade.TEx, blade.UGx];
    y_comb = [blade.y_suction, blade.LEy, blade.y_pressure, blade.TEy, blade.UGy];
    x_ss_comb = flip([blade.UGx, blade.x_suction]);
    y_ss_comb = flip([blade.UGy, blade.y_suction]);
end

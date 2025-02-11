%% GENERATES BLADE GEOMETRY
function [blade, failcode] = pritchard(params, exclusion_factor)
    failcode = "success!!";
    max_iter = 300;

    params.beta_IN = deg2rad(params.beta_IN);
    params.beta_OUT = deg2rad(params.beta_OUT);
    params.ep_IN = deg2rad(params.ep_IN);
    params.ep_OUT = deg2rad(params.ep_OUT);
    params.zeta = deg2rad(params.zeta);   

    o = 2*pi*params.R/params.N_B*cos(params.beta_OUT)-2*params.R_TE;

    factor = 1;
    ttc = params.ttc/100;
    iteration_threshold = 0.001;
    iter_counter = 0;
    while factor > iteration_threshold
        [pts, failcode] = points(params.R, params.R_LE, params.R_TE, params.Cx, params.Ct, params.zeta, params.beta_IN, params.ep_IN, params.beta_OUT, params.ep_OUT, params.N_B, o);
        blade = generate_curves(pts.x_coords, pts.y_coords, pts.betas, params.R_LE, params.R_TE, pts.R_new, o);
        [blade.x_comb, blade.y_comb, blade.x_ss_comb, blade.y_ss_comb] = combiner(blade);     % Combine separate sections into big x and y vectors (good for solidworks exporting later too)
    
        if params.Ct == 0
            params.Ct = pts.Ct;
        end
        [t_max, pos_ss, pos_ps] = max_t(blade);
        chord = sqrt(params.Ct^2 + params.Cx^2);
        % fprintf("TTC = %.3f\n", t_max/chord)
        factor = abs(t_max/chord-ttc);
        % fprintf("factor: %.3f\n", factor)
        if factor > iteration_threshold
            params.Ct = params.Ct*(3+t_max/chord/ttc)/4;
        end

        iter_counter = iter_counter + 1;
        if iter_counter > max_iter
            failcode = "ttc didn't converge";
            break
        end
    end

    % fprintf("Ct: %.3f", params.Ct)
    
    ps_p1x_target = (blade.x(4)+blade.x(5))/2;
    ss_p1x_target = (blade.x(2)+blade.x(3))/2;
    ps_error = ps_p1x_target-blade.ps_p1x;
    ss_error = ss_p1x_target-blade.ss_p1x;
    avg_error = (abs(ps_error) + abs(ss_error))/2;

    ss_del = blade.x(2)-blade.x(3);
    ss_bez_low = blade.x(3) + exclusion_factor * ss_del;
    ss_bez_high = blade.x(2) - exclusion_factor * ss_del;

    ps_del = blade.x(5)-blade.x(4);
    ps_bez_low = blade.x(4) + exclusion_factor * ps_del;
    ps_bez_high = blade.x(5) - exclusion_factor * ps_del;
    ss_is_bad = not(blade.ss_p1x > ss_bez_low && blade.ss_p1x < ss_bez_high);
    ps_is_bad = not(blade.ps_p1x > ps_bez_low && blade.ps_p1x < ps_bez_high);
    if ss_is_bad || ps_is_bad
        failcode = "beziers are cooked";
    end
    % while ss_is_bad || ps_is_bad
    %     params.Ct = params.Ct+avg_error/100;
    %     pts = points(params.R, params.R_LE, params.R_TE, params.Cx, params.Ct, params.zeta, params.beta_IN, params.ep_IN, params.beta_OUT, params.ep_OUT, params.N_B, o);
    %     blade = generate_curves(pts.x_coords, pts.y_coords, pts.betas, params.R_LE, params.R_TE, pts.R_new, o);
    %     [blade.x_comb, blade.y_comb, blade.x_ss_comb, blade.y_ss_comb] = combiner(blade);     % Combine separate sections into big x and y vectors (good for solidworks exporting later too)
    % 
    %     [~, pos_ss, pos_ps] = max_t(blade);
    %     ps_error = ps_p1x_target-blade.ps_p1x;
    %     ss_error = ss_p1x_target-blade.ss_p1x;
    %     prev_error = avg_error;
    %     avg_error = (abs(ps_error) + abs(ss_error))/2;
    % 
    %     % Tells us if we're absolutely dao mei
    %     if avg_error > prev_error
    %         fprintf("really")
    %     end
    %     fprintf("bad\n")
    % 
    %     % Debugging
    %     fprintf("----------------------------------------------\n")
    %     fprintf("Ct: %.3f --> Error: %.3f\n", params.Ct, avg_error)
    %     fprintf("Boundaries -----------------------------------\n")
    %     fprintf("ss: %.3f --> %.3f --> %.3f\n", ss_bez_low, blade.ss_p1x, ss_bez_high)
    %     fprintf("ps: %.3f --> %.3f --> %.3f\n", ps_bez_low, blade.ps_p1x, ps_bez_high)
    %     fprintf("Fixing " + params.name + " beziers...\n")
    % 
    %     ss_is_bad = not(blade.ss_p1x > ss_bez_low && blade.ss_p1x < ss_bez_high);
    %     ps_is_bad = not(blade.ps_p1x > ps_bez_low && blade.ps_p1x < ps_bez_high);
    % end

    blade.x_thicc = [blade.x_pressure(pos_ps), blade.x_ss_comb(pos_ss)];
    blade.y_thicc = [blade.y_pressure(pos_ps), blade.y_ss_comb(pos_ss)];

    % Extra params
    blade.parameters.R = params.R;
    blade.parameters = params;
    blade.parameters.o = o;
    blade.parameters.pitch = 2*pi*params.R/params.N_B;
    blade.parameters.t_max = max_t(blade);
    blade.parameters.t_min = min_t(blade);
    if blade.parameters.t_min < 2.5
        failcode = "blade too thin";
    end
    blade.parameters.zweifel = (4*pi*params.R) / (params.Cx*params.N_B) * sin(params.beta_IN - params.beta_OUT) * cos(params.beta_OUT)/cos(params.beta_IN);
    blade.parameters.blockage_IN = 2*params.R_LE / (blade.parameters.pitch * cos(params.beta_IN))  * 100;
    blade.parameters.blockage_OUT = 2*params.R_TE / (blade.parameters.pitch * cos(params.beta_OUT)) * 100;
    blade.parameters.chord = sqrt(params.Ct^2 + params.Cx^2);
    blade.parameters.calc_ttc = blade.parameters.t_max/blade.parameters.chord;

    blade.parameters.beta_IN = rad2deg(blade.parameters.beta_IN);
    blade.parameters.beta_OUT = rad2deg(blade.parameters.beta_OUT);
    blade.parameters.ep_IN = rad2deg(blade.parameters.ep_IN);
    blade.parameters.ep_OUT = rad2deg(blade.parameters.ep_OUT);
    blade.parameters.zeta = rad2deg(blade.parameters.zeta); 

    % fprintf("zweifel: %.3f\n", blade.parameters.zweifel)
    % fprintf("blade: " + params.name)
    % fprintf("-#-#-#-#-\n\n\n _/~~~\\_ \n  (O_o)  \n \\__|__/ \n    |    \n  _/ \\_  \n_________\n")
end

%% HELPER FUNCTIONS
% Generates bezier curves for suction and pressure side
function [x_suction, y_suction, x_pressure, y_pressure, ss_p1x, ss_p1y, ps_p1x, ps_p1y, k_max_ss, k_max_ps] = generate_bezier(x, y, betas, UGx_length)
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

    P0 = [x(2), y(2)];
    P1 = [ss_p1x, ss_p1y];
    P2 = [x(3), y(3)];
    
    % First derivatives
    dx_dt = 2*(1-t) * (P1(1)-P0(1)) + 2*t*(P2(1)-P1(1));
    dy_dt = 2*(1-t) * (P1(2)-P0(2)) + 2*t*(P2(2)-P1(2));
    
    % Second derivatives
    d2x_dt2 = 2*(P2(1) - 2*P1(1)+P0(1));
    d2y_dt2 = 2*(P2(2) - 2*P1(2)+P0(2));

    k_max_ss = max(abs(abs(dx_dt .* d2y_dt2 - dy_dt .* d2x_dt2) ./ (dx_dt.^2 + dy_dt.^2).^(3/2)));

    % Bezier for pressure side (number of points matches number of combined
    % points for pressure side bezier + UGT arc
    t = linspace(0,1,1000 + UGx_length);
    b0 = (1-t).^2;
    b1 = 2*t.*(1-t);
    b2 = t.^2;
    x_pressure = x(4)*b0 + ps_p1x*b1 + x(5)*b2;
    y_pressure = y(4)*b0 + ps_p1y*b1 + y(5)*b2;

    P0 = [x(4), y(4)];
    P1 = [ps_p1x, ps_p1y];
    P2 = [x(5), y(5)];
    
    % First derivatives
    dx_dt = 2*(1-t) * (P1(1)-P0(1)) + 2*t*(P2(1)-P1(1));
    dy_dt = 2*(1-t) * (P1(2)-P0(2)) + 2*t*(P2(2)-P1(2));
    
    % Second derivatives
    d2x_dt2 = 2*(P2(1) - 2*P1(1)+P0(1));
    d2y_dt2 = 2*(P2(2) - 2*P1(2)+P0(2));

    k_max_ps = max(abs(abs(dx_dt .* d2y_dt2 - dy_dt .* d2x_dt2) ./ (dx_dt.^2 + dy_dt.^2).^(3/2)));
end

% Generates the throat line
function [x_o, y_o] = generate_throat(x, y, betas, o)
    %% Throat Plotting
    o_slope = -1/tan(betas(2));                                 % Gets perpendicular slope
    x_o = [x(2), x(2) + o*cos(atan(o_slope))];     % X-coords for throat plot
    y_o = o_slope.*(x_o-x(2))+y(2);                             % Y-coords for throat plot
end

% Generates the LE, TE, UGT circular arcs
function blade = generate_arcs(x, y, betas, R_LE, R_TE, R_new)
    [blade.LEx, blade.LEy] = arc(x(3) + R_LE * cos(pi/2 - betas(3)), y(3) - R_LE * sin(pi/2 - betas(3)), R_LE, x(3), x(4), y(3), y(4), 1);
    [blade.TEx, blade.TEy] = arc(x(1) + R_TE * cos(pi/2 - betas(1)), y(1) - R_TE * sin(pi/2 - betas(1)), R_TE, x(5), x(1), y(5), y(1), 1);
    [blade.UGx, blade.UGy] = arc(x(6), y(6), R_new, x(1), x(2), y(1), y(2), 10000);
    blade.x = x;
    blade.y = y;
    blade.betas = betas;
end

% Combines all the generation functions
function blade = generate_curves(x, y, betas, R_LE, R_TE, R_new, o)
    blade = generate_arcs(x, y, betas, R_LE, R_TE, R_new);
    [blade.x_suction, blade.y_suction, blade.x_pressure, blade.y_pressure, blade.ss_p1x, blade.ss_p1y, blade.ps_p1x, blade.ps_p1y, blade.k_max_ss, blade.k_max_ps] = generate_bezier(x, y, betas, length(blade.UGx));
    [blade.x_o, blade.y_o] = generate_throat(x, y, betas, o);
end

% Generates Points 1-5 and their betas
function [pts,failcode] = points(R, R_LE, R_TE, Cx, Ct, zeta, beta_IN, ep_IN, beta_OUT, ep_OUT, N_B, o)
    cont = 1;
    max_iter = 100;
    iter_counter = 0;
    failcode = "success!!";
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
        if Ct == 0
            Ct = y2 + (x2-x3)/(beta_2-beta_3) * log(cos(beta_2)/cos(beta_3)) - R_LE * cos(beta_3);
        end
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
        if delta > 0.0001
            ep_OUT = ep_OUT * (y2/yy2);
        else
            cont = 0;
        end

        iter_counter = iter_counter + 1;
        if iter_counter > max_iter
            failcode = "ep_out didn't converge";
            break
        end
    end
    %% Struct Outputs
    pts.betas = [beta_1, beta_2, beta_3, beta_4, beta_5];
    pts.x_coords = [x1, x2, x3, x4, x5, x01];
    pts.y_coords = [y1, y2, y3, y4, y5, y01];
    pts.R_new = R_01;
    pts.Ct = Ct;
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

function [t_min, pos_ss, pos_ps] = min_t(blade)
    thiccnesses = zeros(1, length(blade.x_pressure));
    for i = 1:length(blade.x_pressure)
        thiccnesses(i) = min(vecnorm([blade.x_ss_comb(i); blade.y_ss_comb(i)] - [blade.x_pressure; blade.y_pressure]));
    end
    [t_min, pos_ss] = min(thiccnesses);
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
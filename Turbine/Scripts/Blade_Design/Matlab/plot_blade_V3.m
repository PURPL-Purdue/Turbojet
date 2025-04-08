function plot_blade_V3(blade, plot_throat, plot_t_max, plot_bez_p1, main_linestyle)
    for i = 1:length(blade)
        hold on
        plot3(blade.x_comb,blade.y_comb, ones(1, length(blade.x_comb))*blade.parameters.R, main_linestyle, LineWidth=1)           % Curves!
        plot3(blade.ss_splinex,blade.ss_spliny, ones(1, length(blade.ss_splinex))*blade.parameters.R, '--b', LineWidth=1)           % Curves!
        scatter3(blade.x, blade.y, ones(1, 6)*blade.parameters.R, 20, 'b', "filled")                        % Points!!
        % plot(blade.x_ss_spline_pts, blade.y_ss_spline_pts, 'b--o')
        % fnplt(blade.ss_spline, 'r', 2)

        % papi = spapi(5, blade.x(1:3), blade.y(1:3));
        % fnplt(papi, 'r', 2);
        
        if plot_throat
            line(blade.x_o, blade.y_o, [blade.parameters.R, blade.parameters.R], Color='blue', LineStyle='--')            % Throat!!!
        end
        if plot_t_max
            line(blade.x_thicc, blade.y_thicc, [blade.parameters.R, blade.parameters.R], Color='green', LineStyle='--')   % T-Max!!!!
        end
        if plot_bez_p1
            line([blade.x(3), blade.ss_p1x], [blade.y(3), blade.ss_p1y], [blade.parameters.R, blade.parameters.R], Color='red', LineStyle='--')
            line([blade.ss_p1x, blade.x(2)], [blade.ss_p1y, blade.y(2)], [blade.parameters.R, blade.parameters.R], Color='red', LineStyle='--')
            line([blade.x(4), blade.ps_p1x], [blade.y(4), blade.ps_p1y], [blade.parameters.R, blade.parameters.R], Color='red', LineStyle='--')
            line([blade.ps_p1x, blade.x(5)], [blade.ps_p1y, blade.y(5)], [blade.parameters.R, blade.parameters.R], Color='red', LineStyle='--')
            plot3(blade.ss_p1x, blade.ss_p1y, blade.parameters.R, '.b', MarkerSize=5)
            plot3(blade.ps_p1x, blade.ps_p1y, blade.parameters.R, '.b', MarkerSize=5)
        end

        % plot_verticals(hub, mid, tip, blade_height, R);                                       % Verticals!!!!!
    end
end

function plot_verticals(bladeset_hub, bladeset_mid, bladeset_tip, blade_height, R)
    loft(bladeset_hub, bladeset_mid, blade_height, R, 1)
    loft(bladeset_mid, bladeset_tip, blade_height, R, 2)
end

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

function spliny_the_elder(blade)
    %% Splineless idiots

    % Tangent Line Equations at Points 2-5
    syms x_perp
    % Slopes
    p1_slope = tan(blade.betas(1));
    p2_slope = tan(blade.betas(2));
    p3_slope = tan(blade.betas(3));
    % Tangent Line Equations
    x1_tan_sym = p1_slope*(x_perp-blade.x(1))+blade.y(1);
    x2_tan_sym = p2_slope*(x_perp-blade.x(2))+blade.y(2);
    x3_tan_sym = p3_slope*(x_perp-blade.x(3))+blade.y(3);

    % Intersection of tangent lines to find P1's
    ss_12x = double(vpasolve(x1_tan_sym == x2_tan_sym, x_perp));
    ss_12y = double(subs(x1_tan_sym, x_perp, ss_12x));

    ss_23x = double(vpasolve(x2_tan_sym == x3_tan_sym, x_perp));
    ss_23y = double(subs(x2_tan_sym, x_perp, ss_23x));

    handle_xlength = 1;
    left_handle_x = blade.x(2)-handle_xlength;
    left_handle_y = double(subs(x2_tan_sym, x_perp, left_handle_x));
    right_handle_x = blade.x(2)+handle_xlength;
    right_handle_y = double(subs(x2_tan_sym, x_perp, right_handle_x));

    spline_x = [blade.x(3), ss_23x,     left_handle_x         blade.x(2), blade.x(2), blade.x(2), blade.x(2),    right_handle_x         ss_12x, blade.x(1)]
    spline_y = [blade.y(3), ss_23y,     left_handle_y         blade.y(2), blade.y(2), blade.y(2), blade.y(2),    right_handle_y         ss_12y, blade.y(1)];
    plot(spline_x, spline_y, 'b--o')

    ctrlp = [spline_x; spline_y];
    knots = aptknt(spline_x, 5);
    spline = spmak(knots, ctrlp);

    plot(spline_x, spline_y, 'b--o')
    fnplt(spline, 'r', 2);
end
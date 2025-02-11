function plot_blade_V3(blade, plot_throat, plot_t_max, plot_bez_p1, main_linestyle)
    for i = 1:length(blade)
        plot3(blade.x_comb,blade.y_comb, ones(1, length(blade.x_comb))*blade.parameters.R, main_linestyle)           % Curves!
        scatter3(blade.x, blade.y, ones(1, 6)*blade.parameters.R, 4, 'b', "filled")                        % Points!!
        
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
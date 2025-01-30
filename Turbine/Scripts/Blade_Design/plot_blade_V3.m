function plot_blade_V3(blade)
    for i = 1:length(blade)
        plot3(blade.x_comb,blade.y_comb, ones(1, length(blade.x_comb))*blade.R, '-k')       % Curves!
        scatter3(blade.x, blade.y, ones(1, 6)*blade.R, 4, 'b', "filled")                    % Points!!
        plot3(blade.x_o, blade.y_o, ones(1, length(blade.x_o))*blade.R, '--r')              % Throat!!!
        plot3(blade.x_thicc, blade.y_thicc, (blade.R)*ones(1,20), '--g')                    % T-Max!!!!

        % plot_verticals(hub, mid, tip, blade_height, R);                                   % Verticals!!!!!
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

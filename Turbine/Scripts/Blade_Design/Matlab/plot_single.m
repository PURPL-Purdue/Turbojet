% Modifies the output matricies as required for aesthetics, then plots
function plot_single(blade, plot_throat, plot_t_max, plot_bez_p1, LE_align, num_blade, profiles_to_plot, pitch_align, type)
    % Pitch calculation for multi-blade plotting
    pitch = 2*pi*blade(pitch_align).parameters.R/blade(pitch_align).parameters.N_B;

    % Finding mirroring line
    y_flip = blade(2).parameters.Ct;

    if type == "stator"
        % Flipping Stator upside down and translating back down
        for i = 1:3
            blade(i).y_comb  = 2*y_flip - blade(i).y_comb  - blade(2).parameters.Ct;
            blade(i).y       = 2*y_flip - blade(i).y       - blade(2).parameters.Ct;
            blade(i).y_thicc = 2*y_flip - blade(i).y_thicc - blade(2).parameters.Ct;
            blade(i).y_o     = 2*y_flip - blade(i).y_o     - blade(2).parameters.Ct;
            blade(i).ss_p1y  = 2*y_flip - blade(i).ss_p1y  - blade(2).parameters.Ct;
            blade(i).ps_p1y  = 2*y_flip - blade(i).ps_p1y  - blade(2).parameters.Ct;
        end

        % Lines up leading edges
        if LE_align
            for i = 1:3
                blade(i).y_comb  = blade(i).y_comb  - (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).y       = blade(i).y       - (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).y_thicc = blade(i).y_thicc - (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).y_o     = blade(i).y_o     - (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).ss_p1y  = blade(i).ss_p1y  - (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).ps_p1y  = blade(i).ps_p1y  - (blade(2).parameters.Ct-blade(i).parameters.Ct);
            end
        end
    else
        % Translating rotor down
        for i = 1:3
            blade(i).y_comb  = blade(i).y_comb  - blade(2).parameters.Ct;
            blade(i).y       = blade(i).y       - blade(2).parameters.Ct;
            blade(i).y_thicc = blade(i).y_thicc - blade(2).parameters.Ct;
            blade(i).y_o     = blade(i).y_o     - blade(2).parameters.Ct;
            blade(i).ss_p1y  = blade(i).ss_p1y  - blade(2).parameters.Ct;
            blade(i).ps_p1y  = blade(i).ps_p1y  - blade(2).parameters.Ct;
        end
        % Lines up leading edges
        if LE_align
            for i = 1:3
                blade(i).y_comb  = blade(i).y_comb  + (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).y       = blade(i).y       + (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).y_thicc = blade(i).y_thicc + (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).y_o     = blade(i).y_o     + (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).ss_p1y  = blade(i).ss_p1y  + (blade(2).parameters.Ct-blade(i).parameters.Ct);
                blade(i).ps_p1y  = blade(i).ps_p1y  + (blade(2).parameters.Ct-blade(i).parameters.Ct);
            end
        end
    end
    % Shifting upwards for multiple blade display

    % PLOTTING
    hold on
    [x_low, x_high, y_diff] = scale_graph(blade);
    min_ys = ones(1,length(profiles_to_plot)*(num_blade));
    max_ys = ones(1,length(profiles_to_plot)*(num_blade));

    counter = 1;
    for i = profiles_to_plot            
        plot_blade_V3(blade(i), plot_throat, plot_t_max, plot_bez_p1, '-k')
        [min_ys(counter), max_ys(counter)] = maxmin_y(blade(i));
        counter = counter + 1;
    end
    for j = 1:num_blade-1
        for i = profiles_to_plot
            blade(i).y_comb  = blade(i).y_comb  + pitch;
            blade(i).y       = blade(i).y       + pitch;
            blade(i).y_thicc = blade(i).y_thicc + pitch;
            blade(i).y_o     = blade(i).y_o     + pitch;
            blade(i).ss_p1y  = blade(i).ss_p1y  + pitch;
            blade(i).ps_p1y  = blade(i).ps_p1y  + pitch;
    
            plot_blade_V3(blade(i), plot_throat, plot_t_max, plot_bez_p1, '-r')
            [min_ys(counter), max_ys(counter)] = maxmin_y(blade(i));
            counter = counter + 1;
        end
    end
    y_low = min(min_ys)-2;
    y_high = max(max_ys)+2;

    xlim([x_low,x_high]);
    ylim([y_low,y_high]);
    pbaspect([x_high-x_low,y_high-y_low,blade(2).parameters.blade_height]);
    
end

% Autoscales the plot
function [x_low, x_high, y_diff] = scale_graph(blade)
    x_total = [blade(1).x_comb, blade(2).x_comb, blade(3).x_comb];
    x_low = min(x_total) - 2;
    x_high = max(x_total) + 2;
    y_total = [blade(1).y_comb, blade(2).y_comb, blade(3).y_comb];
    y_low = min(y_total);
    y_high = max(y_total);
    y_diff = y_high - y_low;
end

function [y_low, y_high] = maxmin_y(blade)
    y_low = min(blade.y_comb);
    y_high = max(blade.y_comb);
end
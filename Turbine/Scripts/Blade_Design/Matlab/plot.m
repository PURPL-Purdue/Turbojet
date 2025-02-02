% Modifies the output matricies as required for aesthetics, then plots
function plot(rotor_blade, stator_blade, plot_throat, plot_t_max, plot_bez_p1, LE_align, num_stator, num_rotor, profiles_to_plot, pitch_align)
    
    % Flipping Stator upside down and translating back down
    y_flip = (rotor_blade(2).parameters.Ct + stator_blade(2).parameters.Ct)/2;
    % line([-2,20],[y_flip, y_flip], [stator_blade(2).R, stator_blade(2).R])
    for i = 1:3
        stator_blade(i).y_comb  = 2*y_flip - stator_blade(i).y_comb  - stator_blade(2).parameters.Ct;
        stator_blade(i).y       = 2*y_flip - stator_blade(i).y       - stator_blade(2).parameters.Ct;
        stator_blade(i).y_thicc = 2*y_flip - stator_blade(i).y_thicc - stator_blade(2).parameters.Ct;
        stator_blade(i).y_o     = 2*y_flip - stator_blade(i).y_o     - stator_blade(2).parameters.Ct;
        stator_blade(i).ss_p1y  = 2*y_flip - stator_blade(i).ss_p1y  - stator_blade(2).parameters.Ct;
        stator_blade(i).ps_p1y  = 2*y_flip - stator_blade(i).ps_p1y  - stator_blade(2).parameters.Ct;
    end
    
    % Lines up leading edges
    if LE_align
        for i = 1:3
            stator_blade(i).y_comb  = stator_blade(i).y_comb  - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
            stator_blade(i).y       = stator_blade(i).y       - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
            stator_blade(i).y_thicc = stator_blade(i).y_thicc - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
            stator_blade(i).y_o     = stator_blade(i).y_o     - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
            stator_blade(i).ss_p1y  = stator_blade(i).ss_p1y  - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
            stator_blade(i).ps_p1y  = stator_blade(i).ps_p1y  - (stator_blade(2).parameters.Ct-stator_blade(i).parameters.Ct);
        end
    end

    % Must export at this stage in order to capture flipped stator geometry
    export(rotor_blade, stator_blade)

    % Lines up leading edges
    if LE_align
        for i = 1:3
            rotor_blade(i).y_comb  = rotor_blade(i).y_comb  + (rotor_blade(2).parameters.Ct-rotor_blade(i).parameters.Ct);
            rotor_blade(i).y       = rotor_blade(i).y       + (rotor_blade(2).parameters.Ct-rotor_blade(i).parameters.Ct);
            rotor_blade(i).y_thicc = rotor_blade(i).y_thicc + (rotor_blade(2).parameters.Ct-rotor_blade(i).parameters.Ct);
            rotor_blade(i).y_o     = rotor_blade(i).y_o     + (rotor_blade(2).parameters.Ct-rotor_blade(i).parameters.Ct);
            rotor_blade(i).ss_p1y  = rotor_blade(i).ss_p1y  + (rotor_blade(2).parameters.Ct-rotor_blade(i).parameters.Ct);
            rotor_blade(i).ps_p1y  = rotor_blade(i).ps_p1y  + (rotor_blade(2).parameters.Ct-rotor_blade(i).parameters.Ct);
        end
    end

    % Shifting rotor to the right
    x_offset = 1.2*stator_blade(2).parameters.Cx;
    for i = 1:3
        rotor_blade(i).x_comb  = rotor_blade(i).x_comb  + x_offset;
        rotor_blade(i).x       = rotor_blade(i).x       + x_offset;
        rotor_blade(i).x_o     = rotor_blade(i).x_o     + x_offset;
        rotor_blade(i).x_thicc = rotor_blade(i).x_thicc + x_offset;
        rotor_blade(i).ss_p1x  = rotor_blade(i).ss_p1x  + x_offset;
        rotor_blade(i).ps_p1x  = rotor_blade(i).ps_p1x  + x_offset;
    end

    % Shifting upwards for multiple blade display
    stator_pitch = 2*pi*stator_blade(pitch_align).parameters.R/stator_blade(pitch_align).parameters.N_B;
    rotor_pitch = 2*pi*rotor_blade(pitch_align).parameters.R/rotor_blade(pitch_align).parameters.N_B;

    % PLOTTING
    figure(1)
    hold on
    [x_low, x_high, y_diff] = scale_graph(stator_blade, rotor_blade);
    min_ys = ones(1,length(profiles_to_plot)*(num_stator+num_stator));
    max_ys = ones(1,length(profiles_to_plot)*(num_stator+num_stator));

    counter = 1;
    for i = profiles_to_plot            
        plot_blade_V3(stator_blade(i), plot_throat, plot_t_max, plot_bez_p1)
        [min_ys(counter), max_ys(counter)] = maxmin_y(stator_blade(i));
        counter = counter + 1;
        plot_blade_V3(rotor_blade(i), plot_throat, plot_t_max, plot_bez_p1)
        [min_ys(counter), max_ys(counter)] = maxmin_y(rotor_blade(i));
        counter = counter + 1;
    end
    for j = 1:num_stator-1
        for i = profiles_to_plot
            stator_blade(i).y_comb  = stator_blade(i).y_comb  + stator_pitch;
            stator_blade(i).y       = stator_blade(i).y       + stator_pitch;
            stator_blade(i).y_thicc = stator_blade(i).y_thicc + stator_pitch;
            stator_blade(i).y_o     = stator_blade(i).y_o     + stator_pitch;
            stator_blade(i).ss_p1y  = stator_blade(i).ss_p1y  + stator_pitch;
            stator_blade(i).ps_p1y  = stator_blade(i).ps_p1y  + stator_pitch;
            
            plot_blade_V3(stator_blade(i), plot_throat, plot_t_max, plot_bez_p1)
            [min_ys(counter), max_ys(counter)] = maxmin_y(stator_blade(i));
            counter = counter + 1;
        end
    end
    for j = 1:num_rotor-1
        for i = profiles_to_plot
            rotor_blade(i).y_comb  = rotor_blade(i).y_comb  + rotor_pitch;
            rotor_blade(i).y       = rotor_blade(i).y       + rotor_pitch;
            rotor_blade(i).y_thicc = rotor_blade(i).y_thicc + rotor_pitch;
            rotor_blade(i).y_o     = rotor_blade(i).y_o     + rotor_pitch;
            rotor_blade(i).ss_p1y  = rotor_blade(i).ss_p1y  + rotor_pitch;
            rotor_blade(i).ps_p1y  = rotor_blade(i).ps_p1y  + rotor_pitch;
    
            plot_blade_V3(rotor_blade(i), plot_throat, plot_t_max, plot_bez_p1)
            [min_ys(counter), max_ys(counter)] = maxmin_y(rotor_blade(i));
            counter = counter + 1;
        end
    end
    y_low = min(min_ys)-2;
    y_high = max(max_ys)+2;

    xlim([x_low,x_high]);
    ylim([y_low,y_high]);
    pbaspect([x_high-x_low,y_high-y_low,stator_blade(2).parameters.blade_height]);
    
    triangle_y = y_low + 0.25*(y_diff);
    % plot_triangles(rotor_blade, stator_blade, triangle_vectors)
end

function plot_triangles(rotor_blade, stator_blade, triangle_vectors)
    vel_triangle(V)
end

% Preps the XY matrix and makes it a format that Solidworks likes
function xyz = export_prep(blade, R, Ct)
    z = ones(1,length(blade.x_comb))*R;
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

% Writes the prepped final XY matricies to .txt files
function export(rotor_blade, stator_blade)
    % EXPORTING
    folder = fileparts(fileparts(mfilename('fullpath'))); 

    writematrix(export_prep(stator_blade(1), stator_blade(1).parameters.R, stator_blade(1).parameters.Ct), fullfile(folder, 'Curve_Files\stator_hub.txt'), 'Delimiter', 'tab');
    writematrix(export_prep(stator_blade(2), stator_blade(2).parameters.R, stator_blade(2).parameters.Ct), fullfile(folder, 'Curve_Files\stator_mid.txt'), 'Delimiter', 'tab');
    writematrix(export_prep(stator_blade(3), stator_blade(3).parameters.R, stator_blade(3).parameters.Ct), fullfile(folder, 'Curve_Files\stator_tip.txt'), 'Delimiter', 'tab');

    writematrix(export_prep(rotor_blade(1),  rotor_blade(1).parameters.R,  rotor_blade(1).parameters.Ct),  fullfile(folder, 'Curve_Files\rotor_hub.txt'),  'Delimiter', 'tab');
    writematrix(export_prep(rotor_blade(2),  rotor_blade(2).parameters.R,  rotor_blade(2).parameters.Ct),  fullfile(folder, 'Curve_Files\rotor_mid.txt'),  'Delimiter', 'tab');
    writematrix(export_prep(rotor_blade(3),  rotor_blade(3).parameters.R,  rotor_blade(3).parameters.Ct),  fullfile(folder, 'Curve_Files\rotor_tip.txt'),  'Delimiter', 'tab');
end

% Autoscales the plot
function [x_low, x_high, y_diff] = scale_graph(stator_blade, rotor_blade)
    x_low = min([stator_blade(1).x_comb, stator_blade(2).x_comb, stator_blade(3).x_comb]) - 2;
    x_high = max([rotor_blade(1).x_comb, rotor_blade(2).x_comb, rotor_blade(3).x_comb]) + 2;
    y_total = [stator_blade(1).y_comb, stator_blade(2).y_comb, stator_blade(3).y_comb, rotor_blade(1).y_comb, rotor_blade(2).y_comb, rotor_blade(3).y_comb];
    y_low = min(y_total);
    y_high = max(y_total);
    y_diff = y_high - y_low;
end

function [y_low, y_high] = maxmin_y(blade)
    y_low = min(blade.y_comb);
    y_high = max(blade.y_comb);
end
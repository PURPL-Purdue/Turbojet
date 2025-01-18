function plot_blade(blade, betas, R, R_LE, R_TE, R_new)    
    % Points 1-5
    scatter(blade.points_x, blade.points_y, 2.5, "filled")
    % Suction Side and Pressure Side Polynomials
    plot_polynomials(blade);
    % Throat Line
    plot_throat(blade);
    % Leading Edge, Trailing Edge, and Unguided Turning Section Circles
    plot_arcs(blade.points_x, blade.points_y, betas, R, R_LE, R_TE, R_new)
end

%% Functions For Polynomials, Throat Line, and Circular Arcs
function plot_polynomials(blade)
    plot(blade.x_suction, blade.y_suction, '-b')
    plot(blade.x_pressure,blade.y_pressure, '-r')
end

function plot_throat(blade)
    plot(blade.x_o, blade.y_o, '-k')
end

function plot_arcs(x, y, betas, R, R_LE, R_TE, R_new)
    hold on
    arc(x(3) + R_LE * cos(pi/2 - betas(3)), y(3) - R_LE * sin(pi/2 - betas(3)), R_LE, x(3), x(4), y(3), y(4))
    arc(x(1) + R_TE * cos(pi/2 - betas(1)), y(1) - R_TE * sin(pi/2 - betas(1)), R_TE, x(5), x(1), y(5), y(1))
    % arc(x(1) + R * cos(pi/2 - betas(1)), y(1) - R * sin(pi/2 - betas(1)), R, x(1), x(2), y(1), y(2))
    arc(x(6), y(6), R_new, x(1), x(2), y(1), y(2))
end

%% Helper Function
function arc(x,y,r, x1, x2, y1, y2)
    start = atan2(y1-y,x1-x);
    stop = atan2(y2-y,x2-x);
    if start > stop
        stop = stop + 2 * pi;
    end

    theta = start:pi/10000:stop;
    xunit = r * cos(theta) + x;
    yunit = r * sin(theta) + y;
    
    plot(xunit, yunit, '-k');
end

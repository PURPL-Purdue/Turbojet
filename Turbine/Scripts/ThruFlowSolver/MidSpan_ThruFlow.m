function results = MidSpan_ThruFlow(reaction, blade_height_mm, r_tip_mm, rpm, mass_flux, ...
    inlet_total_pres, inlet_total_temp, backpressure, stator_loss, rotor_loss, show_plot)
    % One Dimensional Mid-Span Flow Solver
    % This script finds the through-flow solution for the turbine using
    % one-dimensional relations and given losses.
    %
    % Author(s): 
    % - Avidh Bavkar [abavkar@purdue.edu]
    %
    % INPUTS:
    % Reaction: Stage Degree of Reaction (0 = Impulse -> 1: Reaction)
    % Blade Height: [in mm] height of blade
    % Tip Radius:   [in mm] radius of the blade tip (turbine outer radius)
    % RPM:          [in rev/min] shaft speed
    % Mass Flux:    [in kg/s] mass flow rate through turbine
    % Inlet Total Pressure: [in pascals] total pressure at turbine inlet
    % Inlet Total Temperature:  [in Kelvin] total temp. at turbine inlet
    % Backpressure: [in pascals] static backpressure behind turbine
    % Stator Loss:  Enthalpy Loss Coefficient of Stator
    % Rotor Loss:   Enthalpy Loss Coefficient of Rotor
    % Show Plot:    set to true to show convergence plot/triangles.


    % I hate writting out the Isentropic Relations a bunch of times, I made
    % a "library" that has these already done so I don't miss a gamma and
    % break everything:
    addpath FlowTools/OneDim_FlowTools/

    
    %% Quantities:
    bar2pa = 10.^5;
    g = 1.4;
    gm1 = g-1;
    gp1 = g+1;
    R = 287; %m^2/s^2K;
    CP = g.*R./gm1;
    
    % Known Quantities:
    r_tip = r_tip_mm./1000; %meters
    r_hub = (r_tip-blade_height_mm./1000); %meters
    % rpm = 80.*1000;
    % mass_flux = 2.09; %kg/s 
    
    % P01 = 300.*1000; %pascals
    P01 = inlet_total_pres;
    % T01 = 1000; %Kelvin
    T01 = inlet_total_temp;
    % P3 = P01.*ts_press_ratio; % Turbine Backpressure.
    P3 = backpressure;

    % Zweifel coeffs used to suggest blade pitch -> blade count. Higher
    % values indicate higher blade loading (1.27 is considered a "high-lift
    % turbine"). The aerodynamics get more critical as blade-loading is
    % increased, if we keep it at 1 (nominal loading) then should be a
    % pretty standard turbine design.
    stator_zweifel_coeff = 1; % 0.8-1 = nominal -> 1.27 = high lift
    rotor_zweifel_coeff  = 1; % 0.8-1 = nominal -> 1.27 = high lift
    
    % To start the iteration:
    massflux_error = 10;
    iteration = 1;
    
    failcode = "success"; %Track if our turbine is invalid for some reason
    max_iterations = 10000; %Maximum iterations to try to converge (if not we fail)
    %% Mean-Line Solution:
    % We can solve everything perfectly without losses, but when we add the
    % losses, the mass-flow is affected. For example a loss in the rotor
    % when solved leads to less mass-flow through the rotor than we have in
    % the stator, which is obviously violating the conservation of mass. We
    % fix this by decreasing the mach number exiting the rotor a little
    % bit. This is why we have to iterate below:
    mach_loss = 0; % In order to conserve mass when we have losses, we lose a little speed.
    while (abs(massflux_error) > 0.0001 && failcode == "success")

        if iteration >= max_iterations
            if(show_plot)
                disp("ERROR: Failed to converge :(");
            end
            failcode = "massflux_convergence_failed";
        end

        % Find the radii & annular area (A1):
        r_mid = (r_tip+r_hub)./2;
        A1 = pi.*(r_tip.^2-r_hub.^2);

        % From the given reaction fraction, we can determine P2 from P3 and P01
        P2 = P01.*(1-(1-reaction).* (1 - (P3./P01).^(gm1/g))).^(g/gm1);
        
        %validate the algebra (solve for reaction):
        actual_reaction = 1 - (1-IsentropicFlow_Tools.tratio_from_pratio(g, (P2/P01)))./...
                              (1-IsentropicFlow_Tools.tratio_from_pratio(g, (P3./P01)));
    
        if abs(reaction - actual_reaction)>10*eps
            % This shouldnt happen unless you mess with some pressures
            if (show_plot)
                disp("ERROR: reactions wrong");
            end
            failcode = "reaction_error";
        end
        
        % Stator Thermodynamics:
        P02s = P01;
        T02s = T01;
        
        M2s = IsentropicFlow_Tools.mach_FROM_totstat_pratio(g, (P02s/P2));
        T2s = T02s./IsentropicFlow_Tools.tratio_from_pratio(g, (P02s/P2));
        
        % Isentropic loss coefficient (not used, enthalpy instead cuz
        % soderberg epic and better):
        %         v2s = M2s.*sqrt(g.*R.*T2s);
        %         v2 = sqrt(v2s.^2.*(1-stator_loss)); %Account for loss.
        
        %Assume the stator is adiabatic:
        T02 = T02s;
        
        %         T2 = T02 - v2.^2./(2.*CP); 
        T2 = (T2s + stator_loss.*T02s)./(1+stator_loss);
        v2 = sqrt(2.*CP.*(T02-T2));
        P02 = P2.*IsentropicFlow_Tools.pratio_from_tratio(g, (T02./T2));
    
        % Compute A2 from mass-flow:
        rho2 = P2./(R.*T2); %Ideal gas law;
    
        A2 = mass_flux./(rho2.*v2); %Consv. of mass
    
        alpha2 = acos(A2./A1); %Trig.
        
        % Stator Vel. Triangle:
        u2 = (rpm.*2.*pi./60).*r_mid;
        U2 = [0, u2];
        V2 = v2.*[cos(alpha2) , sin(alpha2)];
        W2 = V2-U2;
        w2 = norm(W2);
       
        alpha1 = 0; %assumed axial inflow
        beta2  = -atan(W2(2)./W2(1));
        
        % Rotor Thermodynamics:
        % Transform quantities to relative frame:
        M2R = w2./sqrt(g.*R.*T2);
        P02R = P2.*IsentropicFlow_Tools.totstat_pratio_FROM_mach(g, M2R);
        T02R = T2.*IsentropicFlow_Tools.totstat_tratio_FROM_mach(g, M2R);
        
        P03Rs = P02R;
        T03Rs = T02R;
        
        T3s = T03Rs./IsentropicFlow_Tools.tratio_from_pratio(g, (P03Rs./P3));
    
        M3Rs = IsentropicFlow_Tools.mach_FROM_totstat_tratio(g, (T03Rs./T3s));
        w3s = M3Rs.*sqrt(g.*R.*T3s);
                
        % Assume the rotor is also adiabatic:
        T03R = T03Rs;
    
        T3 = (T3s+rotor_loss.*T02R)./(1+rotor_loss);
    
        P03R = P3.*(T03R/T3).^(g/gm1);
    
        M3R = IsentropicFlow_Tools.mach_FROM_totstat_pratio(g, (P03R./P3))+mach_loss;
        w3 = M3R.*sqrt(g.*R.*T3);
        
        % Area/Flow Turning:
        A3_to_A2 = IsentropicFlow_Tools.aratio_FROM_mach(g, M3R)./...
                   IsentropicFlow_Tools.aratio_FROM_mach(g, M2R);
        
        beta3 = acos(cos(beta2).*(A3_to_A2));
        
        % Rotor Vel. Triangle:
        U3 = U2; %Since the radii do not change
        u3 = u2; %Since the radii do not change
        
        W3 = w3 .* [cos(beta3), -sin(beta3)];
        W3s = w3s .* [cos(beta3), -sin(beta3)];
        V3 = W3 + U3;
        V3s = W3s + U3;
        
        v3 = sqrt(V3(1).^2+V3(2).^2);
        v3s = sqrt(V3s(1).^2+V3s(2).^2);
    
        T03 = T3 + v3.^2./(2.*CP);
        T03s = T3s + v3s.^2./(2.*CP);

        P03 = P3.*IsentropicFlow_Tools.pratio_from_tratio(g, T03./T3);
        
        % Check the mass-flow convergence:
        annular_area = pi.*(r_tip.^2-r_hub.^2);
        
        rho2 = P2./(R.*T2);
        stator_massflux = rho2.*V2(1).*annular_area; %First component of V2 is axial.
        
        rho3 = P3./(R.*T3);
        rotor_massflux = rho3.*V3(1).*annular_area;  %First component of V3 is axial.
        
        massflux_error = rotor_massflux-stator_massflux;

        % Check to see if we made a supersonic turbine (would be epic tho):
        M2 = v2./sqrt(g.*R.*T2);
        M3 = v3./sqrt(g.*R.*T3);
        if (M2 > 1) || (M2R > 1) || (M3 > 1) || (M3R > 1)
            if(show_plot)
                disp("ERROR: your turbine is supersonic:(");
            end
            failcode = "supersonic_flow_present";
        end

        % Compute some output/performance parameters:
        power = mass_flux .* (u2.*V2(2) - u3.*V3(2));
        tottot_eff = (T01-T03)./(T01-T03s);
        totstat_eff = (T01-T03)./(T01-T3s);

        stator_ptoc_zweifel = ...
            stator_zweifel_coeff./(2.*cos(-alpha1).^2.*(tan(alpha1)-tan(-alpha2)));

        rotor_ptoc_zweifel = ...
            rotor_zweifel_coeff./(2.*(tan(-beta2)+tan(beta3).*(cos(beta3).^2)));
        
        convergence_quantity = ...
            [iteration, stator_massflux, rotor_massflux, P2, P3, T2, T3, rho2, rho3, v2, v3, w2, w3, power, M2R, M3R];
        if exist('convergence_quantities','var') == 1
            convergence_quantities = [convergence_quantities;convergence_quantity];
        else
            convergence_quantities = convergence_quantity;
        end

        if ~isreal(convergence_quantities) && ~isreal(power)
            if(show_plot)
                disp("ERROR: complex data :(")
            end
            failcode = "complex_data_present";
        end

        if abs(totstat_eff) > 1 || abs(tottot_eff) > 1
            if(show_plot)
                disp("ERROR: over 100% efficiency :(")
            end
            failcode = "nonphysical_efficiencies";
        end
        
        mach_loss = mach_loss - massflux_error./10; %Tune this until your massflux is conserved
        iteration = iteration + 1;
    end
    
    if (show_plot)
        % Draw the velocity triangle :
        figure(1);
        clf;
        title("Mid-Span Triangle");
        draw_vel_triangle(V2, W2, "2", "-");
        draw_vel_triangle(V3, W3, "3", "-");
        axis equal;

        % Plot convergence quantities
        f=figure(2);
        f.Position = [100 100 400 800];
        sgtitle("Rotor Mass-Flux Convergence");
        subplot(7,1,1);
        plot(convergence_quantities(:, 1), convergence_quantities(: ,2), '-k');
        hold on;
        plot(convergence_quantities(:, 1), convergence_quantities(: ,3), '-b');
        ylabel("Mass Flux (kg/s)");
        legend("Inlet", "Outlet");
        
        subplot(7,1,2);
        plot(convergence_quantities(:, 1), convergence_quantities(: ,4)./bar2pa, '-k');
        hold on;
        plot(convergence_quantities(:, 1), convergence_quantities(: ,5)./bar2pa, '-b');
        ylabel("Static Pres. (bar)");
        legend("Inlet", "Outlet");
        
        subplot(7,1,3);
        plot(convergence_quantities(:, 1), convergence_quantities(: ,6), '-k');
        hold on;
        plot(convergence_quantities(:, 1), convergence_quantities(: ,7), '-b');
        ylabel("Static Temp. (K)");
        legend("Inlet", "Outlet");
        
        subplot(7,1,4);
        plot(convergence_quantities(:, 1), convergence_quantities(: ,8), '-k');
        hold on;
        plot(convergence_quantities(:, 1), convergence_quantities(: ,9), '-b');
        ylabel("Density (kg/m^3)");
        legend("Inlet", "Outlet");
        
        subplot(7,1,5);
        plot(convergence_quantities(:, 1), convergence_quantities(: ,10), '-k');
        hold on;
        plot(convergence_quantities(:, 1), convergence_quantities(: ,11), '-b');
        plot(convergence_quantities(:, 1), convergence_quantities(: ,12), '--k');
        plot(convergence_quantities(:, 1), convergence_quantities(: ,13), '--b');
        ylabel("Velocity");
        legend("Inlet Abs.", "Outlet Abs.", "Inlet Rel.", "Outlet Rel.");
        
        subplot(7,1,6);
        plot(convergence_quantities(:, 1), convergence_quantities(: ,14)./1000, '-k');
        hold on;
        ylabel("Kilowatts");
        
        subplot(7,1,7);
        plot(convergence_quantities(:, 1), convergence_quantities(: ,15), '-k');
        hold on;
        plot(convergence_quantities(:, 1), convergence_quantities(: ,16), '-b');
        ylabel("Rel. Mach");
%         ylim([0, 1]);
        legend("Inlet", "Outlet");
    end
    
    %% Output:
    results.U2 = U2;
    results.W2 = W2;
    results.V2 = V2;

    results.U3 = U3;
    results.W3 = W3;
    results.V3 = V3;

    results.alpha1_deg = rad2deg(alpha1);
    results.alpha2_deg = rad2deg(alpha2);
    results.beta2_deg  = rad2deg(beta2);
    results.beta3_deg  = rad2deg(beta3);
    
    results.M2  = M2;
    results.M3  = M3;
    results.M2s = M2s;
    results.M2R = M2R;
    results.M3R = M3R;
    results.P01 = P01;
    results.P02 = P02;
    results.P03 = P03;
    results.P02R = P02R;
    results.P02s = P02s;
    results.P03R = P03R;
    results.P03Rs = P03Rs;
    results.P2 = P2;
    results.P3 = P3;

    results.stator_ptoc_zweifel = stator_ptoc_zweifel;
    results.rotor_ptoc_zweifel  = rotor_ptoc_zweifel;

    results.convergence_quantities = convergence_quantities;
    results.power = power;
    results.massflux = (rotor_massflux+stator_massflux)./2;
    results.mach_loss = mach_loss;
    results.rpm = rpm;
    results.T01 = T01;
    results.T02 = T02;
    results.T02R = T02R;
    results.T02s = T02s;
    results.T03R = T03R;
    results.T2 = T2;
    results.T2s = T2s;
    results.T3 = T3;
    results.T3s = T3s;
    results.tteff = tottot_eff;
    results.tseff = totstat_eff;
    results.fail = failcode;
end

%% Functions
function draw_vel_triangle(V, W, station_num, linestyle)
    hold on;
    axis equal;
    xlabel("Axial Velocity");
    ylabel("Tangential Velocity");
    plot_vector([0,0], V, ("V" + station_num), linestyle+'k');
    plot_vector([0,0], W, ("W" + station_num), linestyle+'r');
    plot_vector(W, V, ("U" + station_num), linestyle+'b');
end

function draw_3D_vel_triangle(V,W, radii, span_fractions, linestyles)    
    for i = 1:length(span_fractions)
        r_idx = cast(span_fractions(i).*(length(V)-1)+1, "uint8");
        linestyle = linestyles(i);
        plot_vector3([0,0], V(r_idx, :), radii(r_idx), linestyle+'k');
        plot_vector3([0,0], W(r_idx, :), radii(r_idx), linestyle+'r');
        plot_vector3(W(r_idx, :), V(r_idx, :), radii(r_idx), linestyle+'b');
    end
    scatter3(V(:, 1), V(:, 2), radii(:), '.k');
    scatter3(W(:, 1), W(:, 2), radii(:), '.r');
end

function plot_vector(start_pos, end_pos, title, linspec)
    diff = end_pos-start_pos;
    quiver(start_pos(1), start_pos(2), diff(1), diff(2), 0, linspec);
    text(start_pos(1) + diff(1)./2, start_pos(2) + diff(2)./2, title);
end

function plot_vector3(start_pos, end_pos, r, linspec)
    diff = end_pos-start_pos;
    quiver3(start_pos(1), start_pos(2), r,  diff(1), diff(2), 0, 0, linspec);
%     textscatter3(start_pos(1) + diff(1)./2, start_pos(2) + diff(2)./2, r, title);
end



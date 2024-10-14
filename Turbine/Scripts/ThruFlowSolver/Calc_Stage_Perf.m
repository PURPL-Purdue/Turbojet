function output = Calc_Stage_Perf(deg_of_reaction, blade_height_mm, ...
    backpressure, rotor_height2chord, rotor_thk2chord, stator_height2chord, stator_thk2chord)
    % Calculate Stage Performance:
    % Given some geometrical parameters, calculates the performance of a
    % turbine stage using Soderberg's loss correlations.
    % 
    % Author(s): 
    % - Avidh Bavkar [abavkar@purdue.edu]
    %
    % Degree of Reaction:   Stage Reaction (0: Impulse -> 1: Reaction)
    % Blade Height: [in mm] Height of blade
    % Backpressure: [in pascals] static backpressure behind turbine 
    % Rotor H2C:    Rotor Height:Chord Ratio
    % Rotor T2C:    Rotor Thickness:Chord Ratio
    % Stator H2C:   Stator Height:Chord Ratio
    % Stator T2C:   Stator Thickness:Chord Ratio

    addpath LossCorrelations/
    addpath BladeStress/
    
    %% Stage Information (Constants):
    r_tip = 127./2000; %meters
    r_hub = r_tip-blade_height_mm./1000; %meters
    r_mid = (r_hub+r_tip)./2;

    rpm = 80.*1000; %RPM
    mass_flux = 0.485; %kg/s
    inlet_total_pres = 400.*1000; %pascals
    inlet_total_temp = 1000; %Kelvin
    
    g = 1.4;
    R = 287; %m^2/s^2K;
    
    blade_height = blade_height_mm./1000;
    
    % Get the actual dimensions of the turbine from the ratios:
    rotor_chord = blade_height./rotor_height2chord;
    stator_chord = blade_height./stator_height2chord;

    stator_thk = stator_chord.*stator_thk2chord;
    rotor_thk = rotor_chord.*rotor_thk2chord;
    
    %% Step 1: Run the Thru-Flow without loss to get blade angles for the loss corr.
    lossless_thru_flow = MidSpan_ThruFlow(deg_of_reaction, blade_height_mm, r_tip.*1000, ...
        rpm, mass_flux, inlet_total_pres, inlet_total_temp, backpressure, 0, 0, false);
    
    % From this we can get the blade angles:
    alpha_1 = 0; %Assumed! (Axial Stator Inflow)
    alpha_2 = abs(atan(lossless_thru_flow.V2(2)./lossless_thru_flow.V2(1)));
    alpha_3 = abs(atan(lossless_thru_flow.V3(2)./lossless_thru_flow.V3(1)));
    
    % And the reynolds numbers:
    stator_RE = BladeReynolds( ...
        lossless_thru_flow.P2, ...
        lossless_thru_flow.T2, ...
        R, ...
        norm(lossless_thru_flow.V2), ...
        stator_chord);
    
    rotor_RE = BladeReynolds( ...
        lossless_thru_flow.P3, ...
        lossless_thru_flow.T3, ...
        R, ...
        norm(lossless_thru_flow.V3), ...
        rotor_chord);
    
    stator_AR = (blade_height/stator_chord); %Stator Aspect Ratio (Height:Chord)
    rotor_AR  = (blade_height/rotor_chord);  %Rotor  Aspect Ratio (Height:Chord
    
    stator_TC = (stator_thk/stator_chord); %Stator Thickness:Chord
    rotor_TC  = (rotor_thk/rotor_chord);   %Rotor  Thickness:Chord
    
    %% Step 2: Determine Losses using Solderberg Correlation:
    stator_enthalphy_loss = Solderberg(alpha_1, alpha_2, stator_AR, stator_TC, stator_RE);
    rotor_enthalphy_loss  = Solderberg(alpha_2, alpha_3,  rotor_AR,  rotor_TC,  rotor_RE);
    
    %% Step 3: Run the Mean-Flow Solver again with actual losses:
    thru_flow_sol = MidSpan_ThruFlow(deg_of_reaction, blade_height_mm, r_tip.*1000, ...
        rpm, mass_flux, inlet_total_pres, inlet_total_temp, backpressure, ...
        stator_enthalphy_loss, rotor_enthalphy_loss, false);
    
    fail = thru_flow_sol.fail;

    %% Step 4: Check Blade Stress
    stress_frac = BladeStress(blade_height_mm, rpm, r_hub.*1000);
    if (stress_frac >= 1)
        fail = "blade_yield_stress_exceeded";
    end

    %% Step 5: Output Variables
    output.power = thru_flow_sol.power;
    output.tt_eff = thru_flow_sol.tteff;
    output.ts_eff = thru_flow_sol.tseff;
    output.sol = thru_flow_sol;
    output.stator_loss = stator_enthalphy_loss;
    output.rotor_loss  = rotor_enthalphy_loss;
    output.stress_frac = stress_frac;
    
    output.bladeheight_mm = blade_height_mm;

    output.stator_chord_mm = stator_chord.*1000;
    output.rotor_chord_mm = rotor_chord.*1000;

    output.max_stator_thk_mm = stator_thk.*1000;
    output.max_rotor_thk_mm  = rotor_thk.*1000;

    output.stator_pitch_mm = stator_chord.*thru_flow_sol.stator_ptoc_zweifel.*1000;
    output.rotor_pitch_mm  = rotor_chord.*thru_flow_sol.rotor_ptoc_zweifel.*1000;

    output.stator_blade_count = 2.*pi.*(r_mid.*1000)./output.stator_pitch_mm;
    output.rotor_blade_count = 2.*pi.*(r_mid.*1000)./output.rotor_pitch_mm;

    output.r_hub_mm = r_hub.*1000;
    output.r_tip_mm = r_tip.*1000;

    output.fail = fail;
end
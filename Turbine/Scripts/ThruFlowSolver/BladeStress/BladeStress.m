% Author: Dev Patel
function [stressRatio] = BladeStress(len, rpm, r_shaft)
    %inconel properties
    stress_yield = 800; %MPa
    density = 0.00819; %g/cm^3
    w = rpm * 2 * pi ./ 60; %angular velocity (rad/s)
    
    stress_axial = 1e-6 * density * len * w^2 * (r_shaft + len/2); %MPa

    stressRatio = stress_axial/stress_yield;

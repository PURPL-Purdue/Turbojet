function S = Zweifel(r_m, chord, alpha_1, alpha_2)
    S = chord/2*0.85/((tan(alpha_2)-tan(alpha_1))*cos(alpha_2)^2);
    N = 2*pi*r_m/S;
    N = round(N);
    S = 2*pi*r_m/N;
end
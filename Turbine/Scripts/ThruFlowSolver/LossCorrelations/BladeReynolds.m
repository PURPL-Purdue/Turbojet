function Re = BladeReynolds(pressure, temperature, R, velocity, chord)
    mu_0 = 1.716*10^-5;
    T_0 = 273;
    S_mu = 111;
    mu = mu_0 * (temperature/T_0)^1.5*(T_0+S_mu)/(temperature+S_mu);
    rho = pressure/(R*temperature);
    Re = rho*chord*velocity/mu;
end
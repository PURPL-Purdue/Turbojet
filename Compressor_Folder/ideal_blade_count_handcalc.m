
d2 = 110 / 1000; % mm
U1 = 203.2; % m/s
U2 = 460.8; % m/s
CU1 = 44.7; % m/s
CU2 = 282.4; % m/s
Lb = 51.81 / 1000; % mm
W1 = 200.4; % m/s
W2 = 195.8; % m/s

Ib = (CU2 / U2) - ((U1 * CU1)/(U2^2)); % (eq 4-2)

max_allow_delta_W = 0.45 * (W1 + W2); % (eq 6-22)
delta_W = max_allow_delta_W;


while delta_W > 0
    D_eq = (W1 + W2 + delta_W) / (2 * W2); % diffusion factor (eq 4-41)

    if D_eq < 2
        z = (2 * pi * d2 * U2 * Ib) / (delta_W * Lb); % (eq 4-42)
        break
    else
        delta_W = delta_W - 0.01;
    end
end


fprintf("\nIdeal Blade Count: %0.3f \n", z)
fprintf("Actual blade count needs to be below or equal to ideal blade count \n\n")
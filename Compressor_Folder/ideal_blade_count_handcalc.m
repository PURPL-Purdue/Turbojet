d2 = 110 / 1000; % mm to m
d1 = 72 / 1000; % mm to m
U1 = 203.2; % m/s
U2 = 460.8; % m/s
CU1 = 44.7; % m/s
CU2 = 282.4; % m/s
Lb = 51.81 / 1000; % mm to m
W1 = 200.4; % m/s
W2 = 195.8; % m/s
beta1 = 33.2; % degrees
beta2 = 20.2; % degrees

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


z_max = floor(z);
z_min = 8;

kz = linspace(6.5, 8, 100);

z2 = kz * ((d2 + d1) / (d2 - d1)) * sind((beta1 + beta2) / 2); % Pfleiderer


v1 = [1 2 3];
v2 = [1 2 3];
z_diff = 13; % diffuser blade count
blade_range = [14, 15, 16, 17];

for x = 1:length(blade_range)
    for a = 1 : length(v1)
        for b = 1 :length(v2)
            m(a, b, x) = abs((v1(a) * blade_range(x)) - (v2(b) * z_diff));
        end
    end
end

m14 = m(:,:,1)';
m15 = m(:,:,2)';
m16 = m(:,:,3)';
m17 = m(:,:,4)';

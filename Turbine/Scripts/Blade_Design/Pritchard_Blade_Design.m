clear;clc;clf;

%% Numbers currently set are the example numbers given in the Pritchard paper, and the output matches their example output.
R = 5.5;            % Airfoil radius            | mm
R_LE = 0.031;       % Leading edge radius       | mm
R_TE = 0.016;       % Trailing edge radius      | mm

Cx = 1.102;         % Axial chord               | mm
Ct = 0.591;         % Tangential chord          | mm
zeta = 6.3;         % Unguided turning angle    | degrees
beta_IN = 35;       % Inlet blade angle         | degrees
ep_IN = 8;          % Inlet half wedge angle    | degrees
beta_OUT = -57;     % Exit blade angle          | degrees
ep_OUT = 3.32;      % Exit half wedge angle     | degrees

N_B = 51;           % Number of blades          | N/A
o = 0.337;          % Throat distance           | mm

beta_IN = deg2rad(beta_IN);
beta_OUT = deg2rad(beta_OUT);
ep_IN = deg2rad(ep_IN);
ep_OUT = deg2rad(ep_OUT);
zeta = deg2rad(zeta);

%% POINT ONE - Suction Surface Trailing Edge Tangency Point
beta_1 = beta_OUT - ep_OUT;
x1 = Cx - R_TE * (1 + sin(beta_1));
y1 = R_TE * cos(beta_1);

%% POINT TWO - Suction Surface Throat Point
beta_2 = beta_OUT - ep_OUT + zeta;
x2 = Cx - R_TE + (o + R_TE) * sin(beta_2);
y2 = 2*pi*R/N_B - (o + R_TE) * cos(beta_2);

%% POINT THREE - Suction Surface Leading Edge Tangency Point
beta_3 = beta_IN + ep_IN;
x3 = R_LE * (1 - sin(beta_3));
y3 = Ct + R_LE * cos(beta_3);

%% POINT FOUR - Pressure Surface Leading Edge Tangency Point
beta_4 = beta_IN - ep_IN;
x4 = R_LE * (1 + sin(beta_4));
y4 = Ct - R_LE * cos(beta_4);

%% POINT FIVE - Pressure Surface Trailing Edge Tangency Point
beta_5 = beta_OUT + ep_OUT;
x5 = Cx - R_TE * (1 - sin(beta_5));
y5 = -R_TE * cos(beta_5);

%% SUCTION SIDE POLYNOMIAL
%% Coefficients
d = (tan(beta_3)+tan(beta_2))/(x3-x2)^2 - 2*(y3-y2)/(x3-x2)^3;
c = (y3-y2)/(x3-x2)^2 - tan(beta_2)/(x3-x2) - d*(x3 + 2*x2);
b = tan(beta_2) - 2*c*x2 - 3*d*x2^2;

%% Polynomial
a = y2 - b*x2 - c*x2^2 - d*x2^3;
x_suction = linspace(x3,x2,100);
y_suction = a + b.*x_suction + c*x_suction.^2 + d*x_suction.^3;

%% PRESSURE SIDE POLYNOMIAL
%% Coefficients
d = (tan(beta_4)+tan(beta_5))/(x4-x5)^2 - 2*(y4-y5)/(x4-x5)^3;
c = (y4-y5)/(x4-x5)^2 - tan(beta_5)/(x4-x5) - d*(x4 + 2*x5);
b = tan(beta_5) - 2*c*x5 - 3*d*x5^2;

%% Polynomial
a = y5 - b*x5 - c*x5^2 - d*x5^3;
x_pressure = linspace(x4,x5,100);
y_pressure = a + b.*x_pressure + c*x_pressure.^2 + d*x_pressure.^3;

%% PLOTTING POINTS
figure(1)
hold on

xlim([-2,15])
ylim([-2,15])
pbaspect([1,1,1])

%% Points 1-5
plot(x1, y1, '.b', MarkerSize=10)
plot(x2, y2, '.r', MarkerSize=10)
plot(x3, y3, '.y', MarkerSize=10)
plot(x4, y4, '.g', MarkerSize=10)
plot(x5, y5, '.k', MarkerSize=10)

%% Suction Side and Pressure Side Polynomials
plot(x_suction,y_suction, '-b')
plot(x_pressure,y_pressure, '-r')

%% Leading Edge, Trailing Edge, and Unguided Turning Section Circles
circle(x3 + R_LE * cos(pi/2 - beta_3), y3 - R_LE * sin(pi/2 - beta_3), R_LE)
circle(x1 + R_TE * cos(pi/2 - beta_1), y1 - R_TE * sin(pi/2 - beta_1), R_TE)
circle(x1 + R * cos(pi/2 - beta_1), y1 - R * sin(pi/2 - beta_1), R)

%% Plots circle
function circle(x,y,r)
    hold on
    theta = 0:pi/100:2*pi;
    xunit = r * cos(theta) + x;
    yunit = r * sin(theta) + y;
    plot(xunit, yunit, '-k');
    hold off
end

%Combustor Sizing Code V3
%PURPL Turbojet
%Author: Josh Hyatt
%10/22/2024
%References: 
% https://ejournal.brin.go.id/ijoa/article/download/38/4817/18421 
% ^^Hole size equation is wrong
% https://www.researchgate.net/publication/303557903_Design_and_analysis_of_annular_combustion_chamber_of_a_low_bypass_turbofan_engine_in_a_jet_trainer_aircraft

clear;

% r = 3/39.37;
m_dot = 0.488361;
T03 = 496.677466;
P03 = 400000;
T04 = 1000;
P04 = 376000;

R = 287; %Gas Constant [Nm/kgK]
syms M
gamma_cold = 1.4;
gamma_hot = 1.33;
T3 = T03*(1+(gamma_cold-1)/2*M^2)^(-1);
T4 = T04*(1+(gamma_hot-1)/2*M^2)^(-1);

% A_ref = (R/2*(m_dot_3*T3^.5/P3)^2*pressure_loss/q_ref*(pressure_loss/P3)^-1)^.5; %Casing Reference Area
A_ref = pi*(3/39.37)^2

M_ref3 = vpasolve(m_dot == (P03*(1+(gamma_cold-1)/2*M^2)^(-gamma_cold/(gamma_cold-1))*M*sqrt(gamma_cold*R*T3)*A_ref)/(R*T03*(1+(gamma_cold-1)/2*M^2)^(-1)),M);
T3 = T03*(1+(gamma_cold-1)/2*M_ref3^2)^(-1)
U_ref3 = M_ref3*sqrt(gamma_cold*R*T3)

speed_of_sound_3 = sqrt(gamma_cold*R*T3)

P3 = P03*(1+(gamma_cold-1)/2*M_ref3^2)^(-gamma_cold/(gamma_cold-1));
P4 = P3 * .94;
% P3 = P03*(1+(gamma_hot-1)/2*M_ref2^2)^(-gamma_hot/(gamma_hot-1));

M_ref4 = vpasolve(m_dot == (P4*M*sqrt(gamma_hot*R*T4)*A_ref)/(R*T04*(1+(gamma_hot-1)/2*M^2)^(-1)),M)
U_ref4 = M_ref4*sqrt(gamma_hot*R*T4);

T3 = T03*(1+(gamma_hot-1)/2*M_ref3^2)^(-1);
T4 = T04*(1+(gamma_cold-1)/2*M_ref4^2)^(-1);

rho_3 = P3/(R*T3);

pressure_loss = (P3-P4);
%Add pressure loss equation???

m_dot_3 = 0.488361;
m_dot_4 = 0.495966;

U_ref = m_dot_3 / (rho_3 * A_ref);
q_ref = rho_3 * U_ref^2 / 2;
M_ref = U_ref / (gamma_cold*R*T3)^.5;

diff_loss = 0;

A_h_eff = A_ref / (pressure_loss/q_ref-diff_loss/q_ref)^.5; %total effective area of the holes

T_max = 1100; %max temp [K]
% PF = (T_max - T04)/(T04-T03) %Pattern Factor

PF = 0.25; %Pattern Factor, suggested range 0.25-0.4

K_opt = 0.66;
A_L = K_opt * A_ref; %Liner Area

%Liner diameters
D_out = sqrt(A_L/pi)*2 %outer liner
D_in = 1.5/39.37 %inner liner (based on turbine)
D_L = (D_out - D_in)/2 %assumes no wall thickness!!!! 


A_ann_out = A_ref - A_L; %Outer Annulus Area

L_pz = 3/4 * D_L %Length of Primary Zone (coefficient suggested 2/3-3/4)

L_sz = 1/2 * D_L %Length of Secondary Zone

PF = 0.25; %Pattern Factor, suggested range 0.25-0.4 (T_max - combustor exit temp)/(combustor exit temp -Combustor inlet temp)

L_dz = D_L * (3.83 - 11.83*PF + 13.4*PF^2) %Dillution Zone Length 

total_length = L_pz + L_sz + L_dz

%Distribution of Airflow

m_dot_pz = .2/2 * m_dot; %primary zone mass flow rate;
m_dot_sz = .3/2 * m_dot; %secondary zone mass flow rate
m_dot_dz  =.5/2 * m_dot; %dillution zone mass flow rate
m_dot_injector = .5 * m_dot; %mass flow rate of air through injector tubes

%Hole Sizes
d_pz = .002; %primary zone hole diameter (m)
d_sz = .0025; %secondary zone hole diameter (m)
d_dz = .003; %dillution zone hole diameter (m)
d_injector = .01; %injector hole diameter (m)

n_pz = 45;
while n_pz > 44
n_pz = 15.25 * m_dot_pz / sqrt(P3*pressure_loss/T3) / d_pz^2;
d_pz = d_pz + .000001;
end

n_sz = 53;
while n_sz > 52
n_sz = 15.25 * m_dot_sz / sqrt(P3*pressure_loss/T3) / d_sz^2;
d_sz = d_sz + .000001;
end

n_dz = 45;
while n_dz > 44
n_dz = 15.25 * m_dot_dz / sqrt(P3*pressure_loss/T3) / d_dz^2;
d_dz = d_dz + .000001;
end

n_injector = 7;
while n_injector > 6
n_injector = 15.25 * m_dot_injector / sqrt(P3*pressure_loss/T3) / d_injector^2;
d_injector = d_injector + .000001;
end

cd = 0.5; %Discharge coefficient (approximated, need jet velocity to solve)

%Actual hole diameters [in]
% n_pz = vpa(n_pz)
d_pz = vpa(d_pz / (cd^.5))*39.37
% 
% n_sz = vpa(n_sz)
d_sz = vpa(d_sz / (cd^.5)) *39.37
% 
% n_dz = vpa(n_dz)
d_dz = vpa(d_dz / (cd^.5)) *39.37

d_injector = vpa(d_injector / (cd^.5)) *39.37



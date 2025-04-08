clc;clear;clf




c = 1/1000 %input('Enter radial clearance in m : ');
s = input('Enter tooth pitch in m: ');
w = input('Enter tooth width in m: ');
R = input('Enter shaft radius in m: ');
N = input('Enter number of teeth : ');
Pi = input('Enter seal inlet pressure in Pa: ');
Pe = input('Enter seal exit pressure in Pa: ');
rho= input('Enter density in kg/m^3 : ');
mu = input('Enter dynamic viscosity in PaS : ');
A = 2*3.142*R*c;
p = ones(1,N+1); p(1) = Pi; p(N+1) = Pe;
mdot = 0.1*A*sqrt((Pi-Pe)*rho);
error=1;
for i=1:10000
 Re = mdot/(3.142*2*R*mu);
 gamma =((1-6.5*(c/s))- 8.638*(c/s) *(w/s))*(Re+((1-6.5*(c/s))- 8.638*(c/s) *(w/s))^(-
1/(2.454*(c/s)+2.258*(c/s)*(w/s)^1.673)))^(2.454*(c/s)+2.258*(c/s)*(w/s)^1.673);
 cd1 =(1.097-0.0029*w/c)*(1+(44.86*w/c)/Re)^(-0.2157)/sqrt(2);
 tdc =cd1*0.925*gamma^.861;
 p(2) = p(1)-(mdot/A)^2/(2*cd1^2*rho);

 for j = 2:N-1
 p(j+1) =p(j)-(mdot/A)^2/(2*tdc^2*rho);
 end
 mdot1 = A* tdc*sqrt(2*(rho*(p(N)-p(N+1))));
 error = abs((mdot1-mdot)/mdot);
 mdot=(mdot1-mdot)*0.1+mdot;
 if(error<0.0001)
 break;
 end
end
display('leakage rate in kg/s is : ');
mdot
display('The pressure distribution is (in Pa) : ');
p 






















%% Code for later hehe
% %% VELOCITY TRIANGLES
% velocity_triangles.V1_mag = 402;
% velocity_triangles.V2_mag = 350;
% velocity_triangles.V3_mag = 375;
% 
% velocity_triangles.RPM = 80000;
% 
% velocity_triangles.U_mag = velocity_triangles.RPM * 2*pi/60 * (T_R - T_blade_height/2)/1000;
% 
% velocity_triangles.a2 = 45;
% velocity_triangles.B2 = -36;
% velocity_triangles.B3 = -45;
% 
% velocity_triangles.W2_mag = 300;
% velocity_triangles.W3_mag = 800;
% 
% % velocity_triangles.max_mag = max(v1_mag, V2_mag, V3_mag, W2_mag, W3_mag, U_mag);




% The cool little guy who tells me my blades are finished
% -#-#-#-#-
% 
%  _/~~~\_ 
%   (O_o)  
%  \__|__/ 
%     |    
%   _/ \_  
% _________

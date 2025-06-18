function [temp,press,rho]=atmosphere4(Hvector,GeometricFlag)
%function [temp,press,rho,Hgeopvector]=atmosphere4(Hvector,GeometricFlag)
% Standard Atmospheric data based on the 1976 NASA Standard Atmoshere.
% Hvector is a vector of altitudes.
% If Hvector is Geometric altitude set GeometricFlag=1.
% If Hvector is Geopotential altitude set GeometricFlag=0.
% Temp, press, and rho are temperature, pressure and density
% output vectors the same size as Hgeomvector.
% Output vector Hgeopvector is a vector of corresponding geopotential altitudes (ft).
% This atmospheric model is good for altitudes up to 295,000 geopotential ft.
% Ref: Intoduction to Flight Test Engineering by Donald T. Ward and Thomas W. Strganac
% index   Lapse rate   Base Temp     Base Geopo Alt        Base Pressure            Base Density
%   i     Ki(degR/ft)  Ti(degR)        Hi(ft)              P, lbf/ft^2           RHO, slug/ft^3
format long g 
D= [1     -.00356616    518.67               0                   2116.22       0.00237691267925741
    2        0          389.97        36089.239          472.675801650081      0.000706115448911997
    3      .00054864    389.97        65616.798          114.343050672041      0.000170813471460564
    4      .00153619    411.57       104986.878          18.1283133205764      2.56600341257735e-05
    5        0          487.17       154199.475          2.31620845720195      2.76975106424479e-06
    6     -.00109728    487.17       170603.675          1.23219156244977      1.47347009326248e-06
    7     -.00219456    454.17       200131.234          0.38030066501701      4.87168173794687e-07
    8        0          325.17       259186.352        0.0215739175227548      3.86714900013768e-08];
	
R=1716.55;	%ft^2/(sec^2degR)
gamma=1.4;
g0=32.17405;	%ft/sec^2
RE=20926476;	% Radius of the Earth, ft
K=D(:,2);	%degR/ft
T=D(:,3);	%degR
H=D(:,4);	%ft
P=D(:,5);	%lbf/ft^2
RHO=D(:,6);	%slug/ft^3

temp=zeros(size(Hvector));
press=zeros(size(Hvector));
rho=zeros(size(Hvector));
Hgeopvector=zeros(size(Hvector));

% Convert from geometric altitude to geopotental altitude, if necessary.
if GeometricFlag
	Hgeopvector=(RE*Hvector)./(RE+Hvector);
	disp('Convert from geometric altitude to geopotential altitude in feet')
else 
   Hgeopvector=Hvector;
   %disp('Input data is geopotential altitude in feet')
end

ih=length(Hgeopvector);
n1=find(Hgeopvector<=H(2));
n2=find(Hgeopvector<=H(3) & Hgeopvector>H(2));
n3=find(Hgeopvector<=H(4) & Hgeopvector>H(3));
n4=find(Hgeopvector<=H(5) & Hgeopvector>H(4));
n5=find(Hgeopvector<=H(6) & Hgeopvector>H(5));
n6=find(Hgeopvector<=H(7) & Hgeopvector>H(6));
n7=find(Hgeopvector<=H(8) & Hgeopvector>H(7));
n8=find(Hgeopvector<=295000 & Hgeopvector>H(8));
icorrect=length(n1)+length(n2)+length(n3)+length(n4)+length(n5)+length(n6)+length(n7)+length(n8);
if icorrect<ih
	disp('One or more altitutes is above the maximum for this atmospheric model')
	icorrect
	ih
end
% Index 1, Troposphere, K1= -.00356616
if length(n1)>0
	i=1;
	h=Hgeopvector(n1);
	TonTi=1+K(i)*(h-H(i))/T(i);
	temp(n1)=TonTi*T(i);
	PonPi=TonTi.^(-g0/(K(i)*R));
	press(n1)=P(i)*PonPi;
	RonRi=TonTi.^(-g0/(K(i)*R)-1);
	rho(n1)=RHO(i)*RonRi;
end

% Index 2,  K2= 0
if length(n2)>0
	i=2;
	h=Hgeopvector(n2);
	temp(n2)=T(i);
	PonPi=exp(-g0*(h-H(i))/(T(i)*R));
	press(n2)=P(i)*PonPi;
	RonRi=PonPi;
	rho(n2)=RHO(i)*RonRi;
end

% Index 3,  K3= .00054864
if length(n3)>0
	i=3;
	h=Hgeopvector(n3);
	TonTi=1+K(i)*(h-H(i))/T(i);
	temp(n3)=TonTi*T(i);
	PonPi=TonTi.^(-g0/(K(i)*R));
	press(n3)=P(i)*PonPi;
	RonRi=TonTi.^(-g0/(K(i)*R)-1);
	rho(n3)=RHO(i)*RonRi;
end

% Index 4,  K4= .00153619
if length(n4)>0
	i=4;
	h=Hgeopvector(n4);
	TonTi=1+K(i)*(h-H(i))/T(i);
	temp(n4)=TonTi*T(i);
	PonPi=TonTi.^(-g0/(K(i)*R));
	press(n4)=P(i)*PonPi;
	RonRi=TonTi.^(-g0/(K(i)*R)-1);
	rho(n4)=RHO(i)*RonRi;
end

% Index 5,  K5= 0
if length(n5)>0
	i=5;
	h=Hgeopvector(n5);
	temp(n5)=T(i);
	PonPi=exp(-g0*(h-H(i))/(T(i)*R));
	press(n5)=P(i)*PonPi;
	RonRi=PonPi;
	rho(n5)=RHO(i)*RonRi;
end

% Index 6,  K6= -.00109728
if length(n6)>0
	i=6;
	h=Hgeopvector(n6);
	TonTi=1+K(i)*(h-H(i))/T(i);
	temp(n6)=TonTi*T(i);
	PonPi=TonTi.^(-g0/(K(i)*R));
	press(n6)=P(i)*PonPi;
	RonRi=TonTi.^(-g0/(K(i)*R)-1);
	rho(n6)=RHO(i)*RonRi;
end

% Index 7,  K7= -.00219456
if length(n7)>0
	i=7;
	h=Hgeopvector(n7);
	TonTi=1+K(i)*(h-H(i))/T(i);
	temp(n7)=TonTi*T(i);
	PonPi=TonTi.^(-g0/(K(i)*R));
	press(n7)=P(i)*PonPi;
	RonRi=TonTi.^(-g0/(K(i)*R)-1);
	rho(n7)=RHO(i)*RonRi;
end

% Index 8,  K8= 0
if length(n8)>0
	i=8;
	h=Hgeopvector(n8);
	temp(n8)=T(i);
	PonPi=exp(-g0*(h-H(i))/(T(i)*R));
	press(n8)=P(i)*PonPi;
	RonRi=PonPi;
	rho(n8)=RHO(i)*RonRi;
end




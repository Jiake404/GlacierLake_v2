function [LW_in] = longwave_in1(rh, air_T)
%GlacierLake function to calculate incoming longwave radiation using 
%equations from Benedek (2014)

e_oc = 0.952;           %emissivity of overcast sky
n = 0.6;                %cloudiness, unitless between 0 and 1
sigma = 5.67*10^-8;     %stefan boltzmann constant [W/m2-K]
b = 0.484;              %constant
m = 8;                  %constant
p = 4;                  %constant

e_sat = 610.8*exp(22.47*(1 - (273.15/air_T))); %for ice

%calculate e_air (reference vapour pressure)
e_air = e_sat*(rh/100); %/100 to go from % to proportion

%calculate e_cs (clear sky emissivity)
e_cs = 0.23 + b*((e_air/air_T)^(1/m));

%calculate LW_in
LW_in = (e_cs*(1 - n^p) + e_oc*(n^p))*sigma*air_T^4;










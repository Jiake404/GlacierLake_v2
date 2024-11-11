function [LW_in] = longwave_in(rh, air_T,  surface_lambda, model_stage)
%GlacierLake function to calculate incoming longwave radiation using 
%equations from Benedek (2014)

% %constants
% C = 8;                 %exchange coefficient for atmospheric stability [W/m2-K]
% L_v = 2.501*10^6;      %latent heat of vaporization [J/kg]
% c_air = 1.005*10^3;    %specific heat capacity of air [J/kg-K]
% p_atm = 90700;         %atmospheric pressure [Pa] (average of 2010 atm pressure during summer melt season)

e_oc = 0.952;           %emissivity of overcast sky
n = 0.6;                %cloudiness, unitless between 0 and 1
sigma = 5.67*10^-8;     %stefan boltzmann constant [W/m2-K]
b = 0.484;              %constant
m = 8;                  %constant
p = 4;                  %constant

%calculate e_sat dependent on surface cell composition
if model_stage == 2 || model_stage == 4 || (model_stage == 1 && surface_lambda == 0)...
        || (model_stage == 3 && surface_lambda == 0) ||...
        (model_stage == 5 && surface_lambda == 0)
    
    e_sat = 610.8*exp(22.47*(1 - (273.15/air_T))); %for ice
else
    e_sat = 610.8*exp(19.858*(1 - (273.15/air_T))); %for water
end

%calculate e_air (reference vapour pressure)
e_air = e_sat*(rh/100); %/100 to go from % to proportion

%calculate e_cs (clear sky emissivity)
e_cs = 0.23 + b*((e_air/air_T)^(1/m));

%calculate LW_in
LW_in = (e_cs*(1 - n^p) + e_oc*(n^p))*sigma*air_T^4;










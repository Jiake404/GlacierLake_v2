function [snow_lambda, snow_ro, snow_c, snow_T, snow_e, snow_k,...
    snow_z, snow_ice, snow_l_ro] = initialise_snow1(snow_mwe,...
    ro_snow_initial, c_ice, k_ice, k_air, ro_water, T_melt, air_T_out)
%GlacierLake subfunction to initialise snow layer

    %set initial values
    snow_lambda = 0; %initially no melt
    snow_ice = 0; %initially no refrozen melt
    snow_ro = ro_snow_initial;
    snow_l_ro = snow_ro;
    snow_c = c_ice; 

    %if air temperature is below freezing, snow falls at
    %temperature of air. If it is above, it falls at freezing point
    if air_T_out < T_melt
        snow_T = air_T_out;
    else
        snow_T = T_melt;
    end
    
    snow_e = snow_T*snow_mwe*snow_c*ro_water; %snow enthalpy
    snow_k = (ro_snow_initial/ro_water)*k_ice + (1 - ro_snow_initial/ro_water)*k_air; %calculate initial conductivity of ice
    snow_z = snow_mwe*(ro_water/ro_snow_initial);  
    
end












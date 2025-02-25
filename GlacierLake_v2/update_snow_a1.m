function [snow_c, snow_e, snow_lambda, snow_ro, snow_k] = update_snow_a1(...
    c_ice, snow_lambda, precip_out, snow_mwe, snow_e, air_T_out, ro_water,...
    T_melt, ro_snow_initial, snow_ro_prev, snow_mwe_prev, k_ice, ro_ice,...
    b_exp, c_water, snow_l_ro, k_snow_ice_max)
%GlacierLake function to update snow at start of stage (as opposed to
%update_snow_b which updates it at the middle-end of the stage)

    %calculate snow layer heat capacity. No need to incorporate snow_ice
    %here as this falls into (1-snow_lambda) as all c_ice
    snow_c = (c_ice*(1 - snow_lambda) + c_water*snow_lambda);

    %add in new snow if present. falls at air_T unless air_T is above freezing

    if precip_out(2) > 0 %& snow_melt_track < 1
        
%        snow_mwe = snow_mwe + precip_out(2);    
            if air_T_out < T_melt
                snow_e = snow_e + air_T_out*ro_water*precip_out(2)*c_ice;
            else
                snow_e = snow_e + T_melt*ro_water*precip_out(2)*c_ice;
            end

        %update lambda to include new snow

        snow_lambda = (snow_lambda.*snow_mwe)./...
            (snow_mwe + precip_out(2));

        %update density to include new snow using weighted average
        %including snow density at previous timestep. 
        snow_ro = (ro_snow_initial*precip_out(2) + snow_ro_prev*snow_mwe_prev)/...
            (snow_mwe_prev+precip_out(2));
        
        %update conductivity after Essery FSM, limit at
        %threshold given by k_snow_ice_max
        snow_k = min(k_ice*(snow_l_ro/ro_ice)^b_exp, k_snow_ice_max); 
    else
        %updates if no snowfall (pass on)
        snow_ro = snow_ro_prev;
        
        %update conductivity using snow_ro_prev after Essery FSM, limit at
        %threshold given by k_snow_ice_max
        snow_k = min(k_ice*(snow_l_ro/ro_ice)^b_exp, k_snow_ice_max);
        
    end

               

end
















function[snow_ro, snow_z, snow_lambda, snow_T, snow_ice, snow_l_ro] = update_snow_b(snow_e, T_melt,...
    snow_mwe, c_ice, ro_water, Lf, snow_ro, ro_melt_max, t_step,...
    tau_ro, ro_cold_max, snow_ice_prev, snow_lambda_prev, snow_l_ro)
%GlacierLake function to update snow at middle-end of stage (as opposed to
%update_snow_a which updates it at the start of the stage)

    %see if snow melt is occuring 
    if snow_e >= T_melt*snow_mwe*c_ice*ro_water

        %calculate snow lambda
        snow_lambda = (snow_e - ro_water*c_ice*snow_mwe*T_melt)./...
            (Lf*ro_water*snow_mwe);

        snow_T = T_melt; %as at melting point

    %if snow is not melting
    else
        %update surface cell temperature
        snow_T = snow_e/(snow_mwe*ro_water*c_ice);
        snow_lambda = 0;

    end
    
    %calculate snow ice layer if water fraction has decreased but is still
    %above zero
    if snow_lambda < snow_lambda_prev && snow_lambda > 0
        
        snow_ice = snow_ice_prev + (snow_lambda_prev - snow_lambda);
        
    %else if no changes keep snow_ice the same as the previous timestep
    elseif snow_lambda > snow_lambda_prev && (snow_lambda + snow_ice_prev) > 1
        snow_ice = 1 - snow_lambda;
    else
        snow_ice = snow_ice_prev;
    end
    
    %when all snow is snow ice treat snow layer as ice layer until it has
    %completely melted
    if snow_ice >= 1
        
        snow_ice = 1;
        snow_z = snow_mwe; %it would be possible to use ro_ice here instead as the snow is a capping layer. But keeping at ro_water to stay consistent with the rest of the model
        snow_ro = ro_water;
        snow_l_ro = ro_water;
        
    else %if snow is not made up entirely of snow_ice
       
        %update snow density 
        if snow_T >= T_melt && snow_ro < ro_melt_max

            %update wet snow after Essery FSM
            snow_ro = ro_melt_max + (snow_ro - ro_melt_max)*...
                exp(-t_step/tau_ro);
            
            %update density to account for relative proportions of snow,
            %water, and snow ice, call snow_l_ro
            snow_l_ro = min(snow_ro*(1 - snow_lambda - snow_ice) + ro_water*(snow_lambda + snow_ice), ro_water);

        elseif snow_T < T_melt && snow_ro < ro_cold_max

            %update dry snow after Essery FSM
            snow_ro = ro_cold_max + (snow_ro - ro_cold_max)*...
                exp(-t_step/tau_ro);    
            
            %update density to account for relative proportions of snow,
            %water, and snow ice, call snow_l_ro
            snow_l_ro = min(snow_ro*(1 - snow_lambda - snow_ice) + ro_water*(snow_lambda + snow_ice), ro_water);
        end
        
        %calculate snow thickness. Take the maximum of snow depth as
        %calculated via snow density and as calculated via snow_mwe (i.e.
        %snow ice)
        snow_z = max([snow_mwe*(1 - snow_lambda - snow_ice)*(ro_water/snow_l_ro),...
            snow_ice*snow_mwe, 0.01]);
        
    end

  
    
    
    
    
    
    
    
    
    
    
    
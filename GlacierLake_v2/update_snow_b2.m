function[snow_ro, snow_z, snow_T, snow_l_ro] = update_snow_b2(snow_lambda, snow_ice,snow_e, T_melt,...
    snow_mwe, c_ice, ro_water, snow_ro, ro_melt_max, t_step, tau_ro, ro_cold_max,  snow_l_ro)
%GlacierLake function to update snow at middle-end of stage (as opposed to
%update_snow_a which updates it at the start of the stage)
    %see if snow melt is occuring 
    if snow_e >= T_melt*snow_mwe*c_ice*ro_water

        snow_T = T_melt; %as at melting point

    %if snow is not melting
    else
        %update surface cell temperature
        snow_T = snow_e/(snow_mwe*ro_water*c_ice);
        snow_T(snow_T<0)=0;

    end
    

    %when all snow is snow ice treat snow layer as ice layer until it has
    %completely melted
    if snow_ice < 1
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

  
    
    
    
    
    
    
    
    
    
    
    
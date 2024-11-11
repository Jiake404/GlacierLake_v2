function [T, lambda, model_stage, snow_T] = temp_lambda_profile1(enthalpy, total_grid_num,...
    grid_profile, T_melt, c_ice, Lf, model_stage, ro_water, c_water, snow_mwe, snow_e, lake_ind)
%GlacierLake subfunction to calculate temperature and water fraction of
%each cell

    %initialise 
    T = zeros(total_grid_num,1);
    lambda = zeros(total_grid_num,1);

    if model_stage == 1
        
        %using ro_water following Benedek (2014)
        if any(enthalpy > c_ice.*ro_water.*T_melt.*grid_profile) %then water present so model stage goes to 3
            model_stage = 3;
            
            return
        else %no full water cells are present
            T = enthalpy./(grid_profile.*ro_water.*c_ice); 
            lambda(:) = 0;
        end
        
    elseif model_stage == 3 || model_stage == 5
        
        %for information on equations see paper. Ice where enthalpy is
        %below latent heat of fusion. lil is 'lambda_ice_logical'.
        lil = enthalpy < c_ice.*ro_water.*T_melt.*grid_profile; 
        lambda(lil) = 0;

        %slush where enthalpy is between that of water and ice. lsl is 
        %'lambda_slush_logical'
        lsl = enthalpy >= c_ice.*ro_water.*T_melt.*grid_profile...
            & enthalpy <= grid_profile.*ro_water.*(T_melt.*c_ice + Lf); 
        lambda(lsl) = (enthalpy(lsl) - ro_water.*c_ice.*grid_profile(lsl)*T_melt)./...
                (Lf*ro_water.*grid_profile(lsl)); %how far along fusion is cell

        %water where enthalpy is above latent heat of fusion. lwl is
        %'lambda_water_logical'
        lwl = enthalpy > grid_profile.*ro_water.*(T_melt.*c_ice + Lf);
        lambda(lwl) = 1;

        %temperature calculations
        T(lil) = enthalpy(lil)./(grid_profile(lil).*ro_water.*c_ice); %calculate enthalpy for ice
        T(lsl) = T_melt;
        T(lwl) = enthalpy(lwl)./(ro_water.*grid_profile(lwl).*c_water) -...
            (1/c_water)*(c_ice*T_melt + Lf - c_water*T_melt);

        if sum(lil + lsl + lwl) > numel(lil)
            error('Indexing error in temp_lambda_profile subfunction. \n \n')
        end
        
        %change model stage as required. Preventing model_stage change when
        %model_stage = 5. Using lid_break_track == 0 as I'm setting this
        %only when it comes from a lid with snow on it 
        if lambda(1) == 0%  && model_stage == 3 %&& lid_break_track == 0;
            if all(lambda == 0)
                model_stage = 1;
            else
                model_stage = 4;
            end
        end
        
        %set to satisfy outputs
        snow_T = 0;
        
    elseif model_stage == 2
        
        %ICE
        %for information on equations see report. Ice where enthalpy is
        %below latent heat of fusion. lil is 'lambda_ice_logical'.
        lil = enthalpy < c_ice.*ro_water.*T_melt.*grid_profile; 
        lambda(lil) = 0;

        %slush cells lambda. Logical array created first. lsl is 
        %'lambda_slush_logical'
        lsl = enthalpy >= c_ice.*ro_water.*T_melt.*grid_profile...
            & enthalpy <= grid_profile.*ro_water.*(T_melt.*c_ice + Lf); 

        lambda(lsl) = (enthalpy(lsl) - ro_water.*c_ice.*grid_profile(lsl)*T_melt)./...
                (Lf*ro_water.*grid_profile(lsl)); %how far along fusion is cell using lsl?

        %water where enthalpy is above latent heat of fusion. lwl is
        %'lambda_water_logical'
        lwl = enthalpy > grid_profile.*ro_water.*(T_melt.*c_ice + Lf);
        lambda(lwl) = 1;

        %temperature calculations. 
        T(lil) = enthalpy(lil)./(grid_profile(lil).*ro_water.*c_ice); %calculate enthalpy for ice
        T(lsl) = T_melt;
        T(lwl) = enthalpy(lwl)./(ro_water.*grid_profile(lwl).*c_water) -...
            (1/c_water)*(c_ice*T_melt + Lf - c_water*T_melt);

        %SNOW
        %calculate temperature of snow dependent on enthalpy
        
        %if all snow is water
        if snow_e > ro_water*snow_mwe*(c_ice*T_melt + Lf)
            
            %snow has melted, so change model stage
            if model_stage == 2; model_stage = 3;
            elseif model_stage == 4; model_stage = 5;
            end
            
            snow_T = T_melt; %hold temp at T_melt
            
        %if enthalpy is between ice and water
        elseif snow_e >= c_ice*ro_water*T_melt*snow_mwe &&...
                snow_e <=snow_mwe*ro_water*(T_melt*c_ice + Lf)
            
            snow_T = T_melt; %hold melting snow layer at T_melt
            
        %all snow is freezing
        else
            snow_T = snow_e/(snow_mwe*ro_water*c_ice);
        end
   
    elseif model_stage == 4 %there is repetition here with stage 2 but makes it easier to track errors and does not effect performance
        
        %ICE
        %for information on equations see report. Ice where enthalpy is
        %below latent heat of fusion. lil is 'lambda_ice_logical'.
        lil = enthalpy < c_ice.*ro_water.*T_melt.*grid_profile; 
        lambda(lil) = 0;

        %slush cells lambda. Logical array created first. lsl is 
        %'lambda_slush_logical'
        lsl = enthalpy >= c_ice.*ro_water.*T_melt.*grid_profile...
            & enthalpy <= grid_profile.*ro_water.*(T_melt.*c_ice + Lf); 

        lambda(lsl) = (enthalpy(lsl) - ro_water.*c_ice.*grid_profile(lsl)*T_melt)./...
                (Lf*ro_water.*grid_profile(lsl)); %how far along fusion is cell using lsl?

        %water where enthalpy is above latent heat of fusion. lwl is
        %'lambda_water_logical'
        lwl = enthalpy > grid_profile.*ro_water.*(T_melt.*c_ice + Lf);
        lambda(lwl) = 1;

        %temperature calculations. 
        T(lil) = enthalpy(lil)./(grid_profile(lil).*ro_water.*c_ice); %calculate enthalpy for ice
        T(lsl) = T_melt;
        T(lwl) = enthalpy(lwl)./(ro_water.*grid_profile(lwl).*c_water) -...
            (1/c_water)*(c_ice*T_melt + Lf - c_water*T_melt);

        %SNOW
        %calculate temperature of snow dependent on enthalpy, but only if
        %snow layer is active
        if snow_mwe > 0
        
            %if all snow is water
            if snow_e > ro_water*snow_mwe*(c_ice*T_melt + Lf)

                %snow has melted, so change model stage
                if model_stage == 2; model_stage = 3;
                elseif model_stage == 4; model_stage = 5;
                end

                snow_T = T_melt; %hold temp at T_melt

            %if enthalpy is between ice and water
            elseif snow_e >= c_ice*ro_water*T_melt*snow_mwe &&...
                    snow_e <=snow_mwe*ro_water*(T_melt*c_ice + Lf)

                snow_T = T_melt; %hold melting snow layer at T_melt

            %all snow is freezing
            else
                snow_T = snow_e/(snow_mwe*ro_water*c_ice);
            end
        else
            snow_T = 0; %to satisfy matlab's hunger for outputs
            
            %does stage need to change before any snow has fallen/
            if enthalpy(1) > c_ice.*ro_water.*T_melt.*grid_profile(1)
                model_stage = 5;
                if lambda(1)>0&&min(lake_ind)==2&&all(diff(lake_ind)==1)
                    model_stage =3;
                end
            end
            
        end
    else
        
        fprintf('This section of temp_lambda_profile sub function has not been completed. \n \n')
    end
        







end
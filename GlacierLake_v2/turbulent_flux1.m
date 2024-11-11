function [enthalpy, dT, lf_lower, lf_upper] = turbulent_flux1(lake_ind, lake_T_av, air_T, t_step, ro_water,...
    c_water, J, T_melt, lake_depth, Lf, enthalpy, grid_profile, T, c_ice, lake_mode, model_stage)
%GlacierLake subfunction to apply turbulent heat flux 

    %model_stage = 5 indicates main lake
    %model_stage = 5.1 indicates surface lake

    %calculate upper and lower turbulent heat fluxes following 
    %Buzzard (2017). *3600*t_step to get flux per time step
    %if lake is open then lf_upper depends on air temperature, if lake is
    %closed then use T_melt
    if model_stage == 3 || model_stage == 5.1 
        lf_upper = sign(lake_T_av - air_T)*ro_water*c_water*J*abs(lake_T_av - air_T)^(4/3)*3600*t_step;
    elseif model_stage == 5 || model_stage == 4
        lf_upper = sign(lake_T_av - T_melt)*ro_water*c_water*J*abs(lake_T_av - T_melt)^(4/3)*3600*t_step;
    end
    lf_lower = sign(lake_T_av - T(max(lake_ind)+1))*ro_water*c_water*J*abs(lake_T_av - T(max(lake_ind)+1))^(4/3)*3600*t_step; %assumes ice-water interface at T_melt
    
    %homogenize lake_T_av and apply turbulent heat fluxes.
    dT = (-lf_upper - lf_lower)/...
        (ro_water*c_water*lake_depth);
    
    %change path depending on lake mode. Do not want to average out lake
    %water if convection is occuring in lake_mode = 2
    if lake_mode == 1
        lake_T_av = lake_T_av + dT;
        T(lake_ind) = lake_T_av;
        
        %recalculate enthalpy 
        enthalpy(lake_ind) = ro_water.*grid_profile(lake_ind).*...
        (c_ice*T_melt + Lf + c_water.*T(lake_ind) - c_water*T_melt);

        %lower turbulent flux goes into underlying ice
        enthalpy(max(lake_ind) + 1) = enthalpy(max(lake_ind) + 1) + lf_lower;
        
        %upper turbulent flux goes into overlying ice if stage 4 or 5, but
        %not if stage is 5.1 meaning subfunction is being used for surface
        %lake
        if model_stage == 4 || model_stage == 5
            if min(lake_ind) - 1 >= 1 %still prevent just in case
                enthalpy(min(lake_ind) - 1) = enthalpy(min(lake_ind) - 1) + lf_upper;
            end
        end
            
    elseif lake_mode == 2
        %dT and lf_lower required, but nothing else      
        
    end

end

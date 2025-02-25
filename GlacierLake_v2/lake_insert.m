function [enthalpy, lambda, T, k, c, SW_prop, grid_profile, tau, Io, track_b, hydro_bucket, total_grid_num,...
    lake_insert_ind, lake_insert_one] = lake_insert(hydro_bucket, ice_grid_z,hydro_T, ro_water, c_ice, T_melt,...
    Lf, c_water, track_b, t_num, tt, enthalpy, lambda, T, k, c, SW_prop, grid_profile,...
    tau, Io,tau_water,k_water,Io_water, total_grid_num,lake_insert_indpre)
%GlacierLake function to initialise hydrograph

    %initialise index of cells to become water
    lake_in_ind = 1:floor(hydro_bucket/ice_grid_z);
    total_grid_num = total_grid_num + length(lake_in_ind);
    %record the amount of overflow added back to the buckets
    hydro_bucket = hydro_bucket - ice_grid_z.*length(lake_in_ind);

    %set enthalpy
    lake_in_e(lake_in_ind) = ro_water.*ice_grid_z.*...
        (c_ice*T_melt + Lf + c_water.*hydro_T - c_water.*T_melt);


    enthalpy_temp = zeros(length(lake_in_ind), t_num);
    lambda_temp = zeros(length(lake_in_ind), t_num);       
    T_temp = zeros(length(lake_in_ind), t_num);           
    k_temp = zeros(length(lake_in_ind), t_num);             
    c_temp = zeros(length(lake_in_ind), t_num);             
    SW_prop_temp = zeros(length(lake_in_ind), t_num);
    grid_profile_temp = zeros(length(lake_in_ind), t_num);
    tau_temp = ones(length(lake_in_ind), t_num)*tau_water; %originally set as tau_ice so this is followed, however not sure if this is necessary
    Io_temp = ones(length(lake_in_ind), t_num);

    % filling
    enthalpy_temp(lake_in_ind,tt) = lake_in_e(lake_in_ind);
    lambda_temp(lake_in_ind,tt) = 1;       
    T_temp(lake_in_ind,tt) = hydro_T;    
    k_temp(lake_in_ind,tt) = k_water;             
    c_temp(lake_in_ind,tt) = c_water;             
    grid_profile_temp(lake_in_ind,tt) = ice_grid_z;
    tau_temp(lake_in_ind,tt) = tau_water; %originally set as tau_ice so this is followed, however not sure if this is necessary
    Io_temp(lake_in_ind,tt) = Io_water;

    % 
    enthalpy = [enthalpy_temp;enthalpy];
    lambda = [lambda_temp;lambda];
    T = [T_temp;T];
    k = [k_temp;k];
    c = [c_temp;c];
    grid_profile = [grid_profile_temp;grid_profile];
    tau = [tau_temp;tau];
    Io = [Io_temp;Io];
    SW_prop = [SW_prop_temp;SW_prop];
    
    if track_b(4) == 0
        lake_insert_ind = length(lake_in_ind);
    else
        lake_insert_ind = lake_insert_indpre + length(lake_in_ind);
    end

    % 
    lake_insert_one = length(lake_in_ind);

    %initialisations complete and hydrograph in use 
    track_b(4) = 1;
    
end
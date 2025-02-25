%copyright 2019 Robert Law, contact rl491@cam.ac.uk or rob3rtlaw@gmail.com
%freely distributed under the terms of the MIT license, including all
%associated subfunctions

%Originally written for an MPhil at the Scott Polar Research Institute, 
%University of Cambridge

%GlacierLake is a MatLab program to model the year-on-year evolution of
%supraglacial lakes as described in Law et al. (2019 or 2020 or never). See
%https://www.spri.cam.ac.uk/people/law/ for updates. GlacierLake is 
%comprised of 5 model stages: 1 = bare ice, 2 = snow on bare ice below 
%threshold, 3 = lake or surface water is present, 4 = lake and lid with
%snow, 5 = lid melt and breakdown

%useful notes:

%time for datasets is referenced in days into AWS dataset. As default the
%model starts on 01/01/xxxx. If altered it may be necessary to amend
%import_snow and import_hydro

%'wubbold' is a flag for any known *small* issues with code which do not
%affect functionality

%possible improvements

%albedo based on snow age after Essery FSM

%=========================================================================
%                               IMPROVEMENTS
%=========================================================================

%GlacierLake_v2 was improved by Jiake Wu from Sun Yat-sen university in
%2024, contact wujk25@mail2.sysu.edu.cn
%See the modification list for details.

%-------------------------------------------------------------------------

clear;clc
tic

%USER INPUTS
AWS = 'SwissCamp';            % Station name
run_start = [2010,01,01];     % Model initiation time
run_end = [2010,12,31];       % Model termination time
albedo_num = 1;               % Albedo-depth parameteration from 0:Luthje et al(2006); 1: Wu et al.(2024 or 2025 or never)
gcnet_num = 1;                % Weather station from 0: PROMICE or 1: GC-NET

%add directories
working_directory = cd;
addpath(strcat(working_directory,'\Data1\RACMO'))
addpath(strcat(working_directory,'\Data1\PROMICE_GC-NET'))
[clip_hour1,clip_hour2,clip_day1,clip_day2,run_day,run_hour] = readtime(AWS,run_start,run_end);

% Parameter setting
[t_step, model_stage, albedo_mult, precip_mult, albedo_add, ice_grid_num, ice_grid_z,...
    deep_ice_grid_num, deep_ice_grid_z, lower_boundary, low_plot_T, hi_plot_T, ice_lim,...
    snow_threshold, slush_lid_num, new_snow_albedo, snow_albedo_min, max_hydro_input,...
    breakup_threshold, lake_turb_threshold, lake_prof_threshold, slush_lid_threshold,...
    slush_mech_threshold, stage_print_hysteresis, refresh_albedo, tau_snow_cold, tau_snow_melt,...
    output_figs, sp_year, AWS_albedo, basal_SW_distribute,...
    ice_albedo, Io_ice, Io_water, Io_slush, tau_ice, tau_water, tau_slush, tau_snow, J, ro_snow_initial,...
    b_exp, ro_melt_max, ro_cold_max, tau_ro, k_snow_ice_max, T_melt, Lf, ro_water, c_ice, c_water,...
    k_water, k_ice, k_air, e_ice, T_ro_max, ro_ice, ies80, t_total, t_num, total_grid_num, grid_profile,...
    lake_grid_num,hydro_T,ice_T_bottom,ice_T_surface] = cons_output(run_day);

% import AWS data
[SW_down_out, LW_down_out, air_T_out, rh_out, hum_out, pressure_out,...
    wind_speed_out, albedo_out, s_data, time_out] = import_AWS1(AWS,clip_hour1,clip_hour2,run_start,gcnet_num);

%apply albedo alterationsalbedo_out(isnan(albedo_out)) = ice_albedo %very few days affected by this.
albedo_out = albedo_out*albedo_mult + albedo_add;
if any(albedo_out > 1)
    fprintf('albedo_mult or albedo_add means albedo values > 1 have occured. \n \n')
    albedo_out(albedo_out > 1) = 1;
end

% import precipitation data
[precip_out,precip_s] = import_snow1(AWS,clip_day1,clip_day2,run_start,run_end,t_num);

% spin_up
[model_stage, ini_condtions] = spin_up(s_data, precip_s, t_step, ice_T_bottom, ice_T_surface,...
    albedo_mult, model_stage, ice_grid_num, ice_grid_z, deep_ice_grid_num, deep_ice_grid_z,...
        lower_boundary, snow_threshold, new_snow_albedo, snow_albedo_min, refresh_albedo,...
        tau_snow_cold, tau_snow_melt, AWS_albedo, ice_albedo, Io_ice, Io_water, Io_slush,...
        ro_snow_initial, b_exp, ro_melt_max, ro_cold_max, tau_ro, k_snow_ice_max,...
        T_melt, Lf, ro_water, c_ice, c_water, k_water, k_ice, k_air, e_ice, ro_ice, sp_year, gcnet_num);

% [hydro_out] = import_runoff(AWS,clip_day1,clip_day2,run_start,run_end,t_num);
[hydro_out,melt_ice,~,Q_mice] = runoff_prod(SW_down_out, LW_down_out, air_T_out, rh_out, hum_out,...
    pressure_out, wind_speed_out, albedo_out, precip_out, run_day, t_step, albedo_mult, model_stage,...
    ini_condtions, ice_grid_num, ice_grid_z, deep_ice_grid_num, deep_ice_grid_z, lower_boundary,...
    snow_threshold, new_snow_albedo, snow_albedo_min, refresh_albedo, tau_snow_cold, tau_snow_melt,...
    ice_albedo, Io_ice, Io_water, Io_slush, ro_snow_initial, b_exp, ro_melt_max,...
    ro_cold_max, tau_ro, k_snow_ice_max, T_melt, Lf, ro_water, c_ice, c_water, k_water, k_ice,...
    k_air, e_ice, ro_ice, gcnet_num);

hydro_out(find(cumsum(hydro_out) >= max_hydro_input, 1, 'first') - 1:end) = 0;


%INITIALISE
%create arrays
enthalpy = zeros(total_grid_num, t_num);        %(J) enthalpy of cells
lambda = zeros(total_grid_num, t_num);          %water content (0 = no water, 1 = all water)
T = zeros(total_grid_num, t_num);               %(K) temperature
k = zeros(total_grid_num, t_num);               %(W/(m.K)) thermal conductivity
c = zeros(total_grid_num, t_num);               %(J/K) heat capacity
SW_prop = zeros(total_grid_num, t_num);         %SW propagation through ice
Io = zeros(total_grid_num, t_num);              %(/m) fraction of shortwave that can propagate beyond surface
tau = ones(total_grid_num, t_num)*tau_ice;      %(/m) bulk shortwave extinction
q = zeros(1, t_num);                            %(J) energy transfer at the surface
lake_depth = zeros(1, t_num);                   %(m) depth of main lake
surface_lake_depth = zeros(1, t_num);           %(m) depth of lake that forms on lid surface
sub_lake_depth = zeros(1, t_num);               %(m) depth of subsurface lake (beneath lid)
lake_albedo = zeros(1, t_num);                  %empty array to record lake albedo change
lid_thick = zeros(1, t_num);                    %(m) thickness of ice lid
lake_T_av = zeros(1, t_num);                    %(K) average (core) lake temperature
lf_upper = zeros(1, t_num);                     %(J) lake flux upper per time step
lf_lower = zeros(1, t_num);                     %(J) lake flux lower per time step
d_enthalpy = zeros(1, t_num);                   %(J) total flux leaving lake core
sur_lake_T_av = zeros(1, t_num);                %(K) average (core) surface lake temperature
snow_z = zeros(1, t_num);                       %(m) overall snow depth including refrozen water
snow_mwe = zeros(1, t_num);                     %(mwe) depth of each layer in mwe
snow_lambda = zeros(1, t_num);                  %lambda of each layer of the snow pack model (0 = no water, 1 = all water)
snow_ice = zeros(1, t_num);                     %of snow pack layer that is not water, how much is snow and how much is snow-ice (0 = all snow, 1 = all snow-ice)
snow_k = zeros(1, t_num);                       %(W/(m.K)) snow thermal conductivity of each snow layer
snow_c = zeros(1, t_num);                       %(J/K) heat capacity of each snow layer
snow_T = zeros(1, t_num);                       %(K) snow temperature
snow_e = zeros(1, t_num);                       %(J) snow enthalpy
snow_ro = zeros(1, t_num);                      %(kg/m^3)density of each snow layer
snow_l_ro = zeros(1, t_num);                    %(kg/m^3)overall density of snow layer, incl water and snow ice
lake_depth_plot = zeros(1,t_num);
lake_insert_ind = zeros(1,t_num);
lake_insert_one = zeros(1,t_num);
hydro_bucket = 0;                               %(m) stores hydro input before threshold is reached
track_a = zeros(1, t_num);                      %tracking array, keep everything in this one and add rows as required.
                                                %row 1 = model stage
track_b = zeros(4,1);                           %second tracking array, this one does not record across whole t_num                                                
                                                %1 = stage 4 snow threshold set
                                                %2 = stage 4 snow initialised
                                                %3 = stage 5 lid breakup
                                                %4 = stage 3 hydrograph input
qm = zeros(1, t_num);
% input initial condtions
ice_T_surface = ini_condtions(1);
ice_T_bottom = ini_condtions(2);
if model_stage == 2
    snow_lambda(1,1) = ini_condtions(3);
    snow_ro(1,1) = ini_condtions(4);
    snow_c(1,1) = ini_condtions(5);
    snow_T(1,1) = ini_condtions(6);
    snow_e(1,1) = ini_condtions(7);
    snow_k(1,1) = ini_condtions(8);
    snow_z(1,1) = ini_condtions(9);
    snow_ice(1,1) = ini_condtions(10);
    snow_l_ro(1,1) = ini_condtions(11);
    snow_mwe(1,1) = ini_condtions(12);
    snow_albedo = ini_condtions(13);
end
%initialise arrays
T(:,1) = linspace(ice_T_surface, ice_T_bottom, total_grid_num)';
enthalpy(:,1) = c_ice.*ro_water.*T(:,1).*grid_profile(:,1);    %following Benedek (2014), all ice to begin with as default
lambda(:,1) = 0;                                                    %begin as ice as default
track_a(1,1) = model_stage;                                           %to avoid 0 error

%PRE LOOP CALCULATIONS
%find when precipitation passes threshold for inclusion
cumulative_snow = cumsum(precip_out(1,:)); %cumulative sum of precipitation
snow_threshold_pass_a = find(cumulative_snow > snow_threshold, 1, 'first');
hydrograph_start = find(hydro_out > 0, 1, 'first'); %timestep to start hydrograph
resize_tt = hydrograph_start;
hydro_i = 1;
hydro_num = 0;

if isempty(hydrograph_start); hydrograph_start = t_num + 1; end

%warnings
if Io_ice > Io_water || Io_slush > Io_water; fprintf('Warning, edit SW_prop to allow for a lower value of Io_water than Io_ice. \n \n'); end

%% MAIN LOOP
fprintf('Main loop commencing. \n \n')

for tt = 1:t_num - 1
    
    % Hydrological input parameters need to be updated during multi-year operation
    if tt > 1
        if time_out(tt) ~= time_out(tt-1)
            track_a(4) = 0;
            hydro_i = 1;
            hydrograph_start = find(hydro_out(find(time_out == time_out(tt),1,'first'):end) > 0, 1, 'first')+find(time_out == time_out(tt),1,'first')-1;
        end
        if track_a(tt)==3&&track_a(tt-1)~=3
            hydro_bucket = 0;
        end
    end
    
    if tt == 1
        snow_mwe(tt) = snow_mwe(tt)+precip_out(2,tt);
    else
        snow_mwe(tt) = snow_mwe(tt-1)+precip_out(2,tt);
    end

    if model_stage == 3 || model_stage == 5
        snow_mwe(tt) = 0;
    elseif model_stage == 1
        if snow_mwe(tt)>= snow_threshold
            model_stage = 2;
            fprintf('Moving from stage %1.0f to %1.0f. Day %3.0f. \n \n', track_a(1,tt), model_stage, tt*t_step/24);
        end
    end

    %STAGE 1, bare ice
    if model_stage == 1
        
        if gcnet_num == 1 %then calculate LW_in seperately
            [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1,tt), model_stage);
        end

        %surface flux 
        [q(tt)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
            wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
            T(1,tt), ice_albedo, AWS_albedo, Io_ice, e_ice,0);
        
        %apply incoming energy flux to enthalpy array
        enthalpy(1,tt) = enthalpy(1,tt) + q(tt)*3600*t_step; %*3600*t_step for amount in that timestep
        
        %calcualte temperature and lambda for cells based on enthalpy
        [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water);
        
        %only continue if model stage still = 1
        if model_stage == 1
            
            %update thermal conductivity and heat capacity to that of ice
            k(:,tt) = k_ice;
            c(:,tt) = c_ice;
            
            %conduction calculations
            [enthalpy(:,tt+1)] = conduct_update_enthalpy_new(k(:,tt), c(:,tt), T(:,tt),...
                grid_profile(:,tt), enthalpy(:,tt), t_step, total_grid_num, ro_water, lower_boundary,...
                model_stage);
            
            %updates
            grid_profile(:,tt + 1) = grid_profile(:,tt);
            T(1,tt + 1) = T(1,tt); %only top cell gets updated for surface flux calculations
            track_a(1,tt + 1) = 1; %model stage = 1
            
        else

            %print that stage is switching if it has been in stage 1 a
            %sufficient number of days
            if all(track_a(1,((stage_print_hysteresis*24)/t_step):tt) == 1)
                fprintf('Moving from stage 1 to %1.0f. Day %3.0f. \n \n', model_stage, tt*t_step/24)
            end

            track_a(1,tt + 1) = 1;
        end
        
    end

    %STAGE 2, snow on ice
    if  model_stage == 2 

        if gcnet_num == 1 %then calculate LW_in seperately
            [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1,tt), model_stage);
        end

        if track_a(1,tt) ~= 2 %if coming from another stage then initialise 
           
            %initialise snow variables
            [snow_lambda(tt), snow_ro(tt), snow_c(tt), snow_T(tt),snow_e(tt), snow_k(tt), snow_z(tt),...
                snow_ice(tt), snow_l_ro(tt)] = initialise_snow1(snow_mwe(tt), ro_snow_initial,...
                c_ice, k_ice, k_air, ro_water, T_melt, air_T_out(tt));
            
            snow_albedo = new_snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end
        
        elseif track_a(1,tt) == 2 %if within stage 2 then update as normal
            
            if tt == 1
                [snow_c(tt), snow_e(tt), snow_lambda(tt), snow_ro(tt),...
                    snow_k(tt)] = update_snow_a1(c_ice, snow_lambda(tt),...
                    precip_out(:,tt), snow_mwe(tt), snow_e(tt), air_T_out(tt), ro_water,...
                    T_melt, ro_snow_initial, snow_ro(tt), snow_mwe(tt), k_ice,...
                    ro_ice, b_exp, c_water, snow_l_ro(tt), k_snow_ice_max);
            else
                %update snow variables
                [snow_c(tt), snow_e(tt), snow_lambda(tt), snow_ro(tt),...
                    snow_k(tt)] = update_snow_a1(c_ice, snow_lambda(tt),...
                    precip_out(:,tt), snow_mwe(tt), snow_e(tt), air_T_out(tt), ro_water,...
                    T_melt, ro_snow_initial, snow_ro(tt - 1), snow_mwe(tt - 1), k_ice,...
                    ro_ice, b_exp, c_water, snow_l_ro(tt), k_snow_ice_max);
            end
            
        end
        
        %surface flux

        [snow_albedo] = snow_albedo_calc(snow_albedo, precip_out(tt),...
            refresh_albedo, t_step, new_snow_albedo, snow_albedo_min, tau_snow_cold,...
            tau_snow_melt, snow_T(tt), T_melt);
        snow_albedo = snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end
        
        [q(tt)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
            wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
            snow_T(tt), snow_albedo, AWS_albedo, Io_ice, e_ice,0); 
        
        %apply surface flux to surface cell
        snow_e(tt) = snow_e(tt) + q(tt)*t_step*3600; %*t_step... to get energy over timestep
        
        %calcualte temperature and lambda for cells based on enthalpy
        [T(:,tt), lambda(:,tt), model_stage, snow_T(tt)] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water, snow_mwe(tt), snow_e(tt));
         if model_stage == 3
             snow_mwe(tt) = 0;
         end
        %only continue if still model stage 2
        if model_stage == 2
         
            %set thermal conductivity and heat capacity of ice cells
            k(:,tt) = lambda(:,tt)*k_water + (1 - lambda(:,tt))*k_ice;
            c(:,tt) = lambda(:,tt)*c_water + (1 - lambda(:,tt))*c_ice;
            
            %conduction calculations
            [enthalpy_temp] = conduct_update_enthalpy_new(k(:,tt), c(:,tt), T(:,tt),...
                grid_profile(:,tt), enthalpy(:,tt), t_step, total_grid_num, ro_water, lower_boundary,...
                model_stage, snow_T(tt), snow_k(tt), snow_c(tt), snow_z(tt), snow_e(tt));
            
            %seperate enthalpy_temp 
            snow_e(tt + 1) = enthalpy_temp(1);
            enthalpy(:,tt + 1) = enthalpy_temp(2:end);
                
            %second round of snow updates
            if tt > 1
            [snow_ro(tt), snow_z(tt), snow_lambda(tt), snow_T(tt), snow_ice(tt), snow_l_ro(tt)] =...
                update_snow_b(snow_e(tt), T_melt, snow_mwe(tt), c_ice, ro_water,...
                Lf, snow_ro(tt), ro_melt_max, t_step, tau_ro, ro_cold_max,...
                snow_ice(tt-1), snow_lambda(tt-1), snow_l_ro(tt));
            else
                [snow_ro(tt), snow_z(tt), snow_lambda(tt), snow_T(tt), snow_ice(tt), snow_l_ro(tt)] =...
                update_snow_b(snow_e(tt), T_melt, snow_mwe(tt), c_ice, ro_water,...
                Lf, snow_ro(tt), ro_melt_max, t_step, tau_ro, ro_cold_max,...
                0, 0, snow_l_ro(tt));
            end
            
            %updates
            grid_profile(:,tt + 1) = grid_profile(:,tt);
            T(1,tt + 1) = T(1,tt); %only top cell gets updated for surface flux calculations
            snow_T(tt + 1) = snow_T(tt); %for surface flux calculations
            snow_e(:,tt + 1) = snow_e(:,tt);
            snow_lambda(:,tt + 1) = snow_lambda(:,tt);
            snow_mwe(:,tt + 1) = snow_mwe(:,tt);
            snow_ro(:,tt + 1) = snow_ro(:,tt);
            snow_l_ro(:,tt + 1) = snow_l_ro(:,tt);
            snow_z(:,tt + 1) = snow_z(:,tt);
            snow_ice(:,tt + 1) = snow_ice(:,tt);
            track_a(1,tt + 1) = 2; %model stage = 2
            
        else
            
            fprintf('Moving from stage 2 to %1.0f. Day %3.0f. \n \n', model_stage, tt*t_step/24)
        end
       
        
    end

    %STAGE 3, lake

    if model_stage == 3 && track_a(1,tt) ~= 3 %if coming from another stage then only do conduction and necessary set up
        
        if gcnet_num == 1 %then calculate LW_in seperately
            [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1,tt), model_stage);
        end

        %recalcualte temperature and lambda for cells based on enthalpy
        [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water);
        
        %calculate lake index for conduction
        [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
            slush_lid_bot] = lake_index_new(lambda(:,tt), ice_lim, total_grid_num, grid_profile(:,tt),...
            tt, t_step, slush_lid_threshold, slush_lid_num);
   

        %record variables
        lake_depth(tt) = bot_lake_depth;
        lake_depth_plot(tt) = bot_lake_depth;
        lake_T_av(tt) = sum(T(lake_ind,tt).*grid_profile(lake_ind,tt))./sum(grid_profile(lake_ind,tt));
        
        %calculate turbulent heat flux if required
        if lake_depth(tt) >= lake_turb_threshold
            lake_mode = 1; %mode 1 = turbulent flux only, no convection
            
            [enthalpy(:,tt)] = turbulent_flux1(lake_ind, lake_T_av(tt), air_T_out(tt),...
                t_step, ro_water, c_water, J, T_melt, lake_depth(tt), Lf,...
                enthalpy(:,tt), grid_profile(:,tt), T(:,tt), c_ice, lake_mode, model_stage);
        end
        
        %recalcualte temperature and lambda for cells based on enthalpy
        [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water);
        
        %set thermal conductivity and heat capacity
        k(:,tt) = lambda(:,tt)*k_water + (1 - lambda(:,tt))*k_ice;
        c(:,tt) = lambda(:,tt)*c_water + (1 - lambda(:,tt))*c_ice;
        
        %conduction calculations
        [enthalpy(:,tt+1)] = conduct_update_enthalpy_new(k(:,tt), c(:,tt), T(:,tt),...
            grid_profile(:,tt), enthalpy(:,tt), t_step, total_grid_num, ro_water, lower_boundary,...
            model_stage);
        
        %updates
        grid_profile(:,tt + 1) = grid_profile(:,tt);
        T(1,tt + 1) = T(1,tt); %only top cell gets updated for surface flux calculations
        track_a(1,tt + 1) = 3; %model stage = 3
        lambda(:,tt + 1) = lambda(:,tt);
        
    elseif model_stage == 3 && track_a(1,tt) == 3 %main stage 3 section (not coming from another stage in same time step)
       
       if tt >= hydrograph_start
            
           if hydro_i == 1

               if hydro_num == 0
                   hydro_bucket = hydro_bucket + sum(hydro_out(1:tt-1));
               else
                   hydro_bucket = hydro_bucket + sum(hydro_out(find(time_out == time_out(tt),1,'first'):tt-1));
               end
               hydro_i = 0;

           end

            %update hydro_bucket
            hydro_bucket = hydro_bucket + hydro_out(tt);
            
            %if hydro_bucket exceeds thickness of a grid cell andhydro_bucket >= ice_grid_z && track_b(4) == 1
            %hydrograph hasn't been used yet, initialise hydrograph
            if hydro_bucket >= ice_grid_z

                % track_b
                [enthalpy, lambda, T, k, c, SW_prop, grid_profile, tau, Io, track_b, hydro_bucket,total_grid_num,...
                    lake_insert_ind(tt),lake_insert_one(tt)]...
                    = lake_insert(hydro_bucket, ice_grid_z,hydro_T, ro_water, c_ice, T_melt,...
                    Lf, c_water, track_b, t_num, tt, enthalpy, lambda, T, k, c, SW_prop, ...
                    grid_profile, tau, Io,tau_water,k_water,Io_water, total_grid_num,lake_insert_ind(tt-1));

            end
        end


        [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water, snow_mwe(tt), snow_e(tt), lake_ind);

        [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
            slush_lid_bot] = lake_index_new(lambda(:,tt), ice_lim, total_grid_num, grid_profile(:,tt),...
            tt, t_step, slush_lid_threshold, slush_lid_num);

         if ~isempty(top_lake_ind)
            lake_depth(tt) = top_lake_depth+bot_lake_depth;
        else
            lake_depth(tt) = bot_lake_depth;
        end

        %calculate lake albedo using lake depth from Benedek (2014)
        [lake_albedo(tt)] = cal_lakealbedo(albedo_num,lake_depth(tt));
        %albedo multiplier 
        lake_albedo(tt) = lake_albedo(tt)*albedo_mult; if lake_albedo(tt) > 1; lake_albedo(tt) = 1; end
        
        if gcnet_num == 1 %then calculate LW_in seperately
            [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1,tt), model_stage);
        end

        %surface flux
        [q(tt),~,~,SWL(tt),~,~,~,q1(tt)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
            wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
            T(1,tt), lake_albedo(tt), 0, Io_water, e_ice,3); %AWS_albedo = 0 as lake used

        %calculate Io and tau profiles
        [Io(:,tt)] = Io_calc(Io_water, Io_ice, Io_slush, total_grid_num, lambda(:,tt));
        [tau(:,tt)] = tau_calc(tau_water, tau_slush, tau_ice, total_grid_num, lambda(:,tt));
        
        %SW propagation through water. 0 is to prevent AWS albedo from
        %being used
        [SW_prop(:,tt)] = SW_propagate1(SW_down_out(tt), grid_profile(:,tt), albedo_out(tt),...
            t_step, Io(:,tt), tau(:,tt), lake_albedo(tt), 0, basal_SW_distribute);
        
        %apply enthalpy changes
        enthalpy(1,tt) = enthalpy(1,tt) + q(tt)*3600*t_step; %*3600*t_step for amount in that timestep
        enthalpy(:,tt) = enthalpy(:,tt) + SW_prop(:,tt);

       [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water, snow_mwe(tt), snow_e(tt), lake_ind);

        % 
        [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
            slush_lid_bot] = lake_index_new(lambda(:,tt), ice_lim, total_grid_num, grid_profile(:,tt),...
            tt, t_step, slush_lid_threshold, slush_lid_num);
        
        if ~isempty(top_lake_ind)
            lake_depth(tt) = top_lake_depth+bot_lake_depth;
%             lake_ind = top_lake_ind;
%             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
        else
            lake_depth(tt) = bot_lake_depth;
%             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
        end

        %otherwise calculate as normal
        lake_T_av(tt) = sum(T(lake_ind,tt).*grid_profile(lake_ind,tt))./sum(grid_profile(lake_ind,tt));
        
        %if lake is above threshold for turbulence calculate flux, calculate fluxes
        %homogenize core temperature if above turbulence threshold, apply
        %turbulence fluxes and update enthalpy
        if lake_depth(tt) >= lake_turb_threshold
            
            %calculate upper and lower turbulent heat fluxes following 
            %Buzzard (2017)
            if min(lake_ind) == 1
                lf_upper(tt) = sign(lake_T_av(tt) - air_T_out(tt))*ro_water*c_water*J*abs(lake_T_av(tt) - air_T_out(tt))^(4/3)*3600*t_step;
            else
                lf_upper(tt) = sign(lake_T_av(tt) - T_melt)*ro_water*c_water*J*abs(lake_T_av(tt) - T_melt)^(4/3)*3600*t_step;
            end
            lf_lower(tt) = sign(lake_T_av(tt) - T(max(lake_ind)+1,tt))*ro_water*c_water*J*abs(lake_T_av(tt) - T(max(lake_ind)+1,tt))^(4/3)*3600*t_step; %assumes ice-water interface at T_melt
            
            %calculate change in temperature of core
            dT = (-lf_upper(tt) - lf_lower(tt))/...
                (ro_water*c_water*lake_depth(tt));

            %apply turbulent fluxes as required
            if min(lake_ind) > 1 %if statement to apply top flux if lid has formed
                enthalpy(min(lake_ind) - 1,tt) = enthalpy(min(lake_ind) - 1,tt) + lf_upper(tt);
            end
            enthalpy(max(lake_ind)+1,tt) = enthalpy(max(lake_ind)+1,tt) + lf_lower(tt);

            [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
                total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
                ro_water, c_water, snow_mwe(tt), snow_e(tt), lake_ind);

            %update core temperature
            lake_T_av(tt) = lake_T_av(tt) + dT;
            T(lake_ind,tt) = lake_T_av(tt);
        end
        
        if lake_depth(tt) >= lake_prof_threshold

            %calculate convection profile and assign to main temperature and 
            %enthalpy profiles
            [T(lake_ind,tt), enthalpy(lake_ind,tt)] = convection1(T(lake_ind,tt),...
                grid_profile(lake_ind,tt), ies80, ro_water, c_ice, T_melt, Lf, c_water);
            
        end
        
        %========================test=================================
        [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
            slush_lid_bot] = lake_index_new(lambda(:,tt), ice_lim, total_grid_num, grid_profile(:,tt),...
            tt, t_step, slush_lid_threshold, slush_lid_num);
        
        if ~isempty(top_lake_ind)
            lake_depth(tt) = top_lake_depth+bot_lake_depth;
%             lake_ind = top_lake_ind;
%             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
        else
            lake_depth(tt) = bot_lake_depth;
%             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
        end
        %========================test=================================

        %set thermal conductivity and heat capacity
        k(:,tt) = lambda(:,tt)*k_water + (1 - lambda(:,tt))*k_ice;
        c(:,tt) = lambda(:,tt)*c_water + (1 - lambda(:,tt))*c_ice;
        
        %conduction calculations
        [enthalpy(:,tt+1)] = conduct_update_enthalpy_new(k(:,tt), c(:,tt), T(:,tt),...
            grid_profile(:,tt), enthalpy(:,tt), t_step, total_grid_num, ro_water, lower_boundary,...
            model_stage, snow_T(tt), snow_k(tt), snow_c(tt), snow_z(tt), snow_e(tt));
        
        %updates
        grid_profile(:,tt + 1) = grid_profile(:,tt);
        T(1,tt + 1) = T(1,tt); %only top cell gets updated for surface flux calculations
        lambda(:,tt + 1) = lambda(:,tt);
        track_a(1,tt + 1) = 3; %model stage = 3 
        lake_insert_ind(1,tt+1) = lake_insert_ind(tt);
        %display stage change if sufficient days have passed 
        if model_stage ~= 3
            stage_1_end = find(track_a(1,:) == 1, 1, 'last');
            if ((tt - stage_1_end)*t_step)/24 >= stage_print_hysteresis
                fprintf('Moving from stage 3 to %1.0f. Day %3.0f. \n \n', model_stage, tt*t_step/24)
            end  
        end
    end 
    
    %STAGE 4, lake with lid
    if model_stage == 4
        if gcnet_num == 1 %then calculate LW_in seperately
            [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1,tt), model_stage);
        end

        %when the moment comes, and everything is in place, unleash the snow
        if snow_mwe(tt)>= snow_threshold && track_b(2) == 0

            %initialise snow variables
            [snow_lambda(tt), snow_ro(tt), snow_c(tt), snow_T(tt),...
                snow_e(tt), snow_k(tt), snow_z(tt), snow_ice(tt), snow_l_ro(tt)] = initialise_snow1(snow_mwe(tt), ...
                ro_snow_initial,c_ice, k_ice, k_air, ro_water, T_melt, air_T_out(tt));
            
            snow_albedo = new_snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end
            track_b(2) = 1; %update tracker

        elseif snow_mwe(tt)>= snow_threshold && track_b(2) == 1 %if the snow has been unleashed, proceed to handle it, but proceed with care #bangers (https://olivercoates.bandcamp.com/album/shelleys-on-zenn-la)

            %update snow variables
            [snow_c(tt), snow_e(tt), snow_lambda(tt), snow_ro(tt),...
                snow_k(tt)] = update_snow_a1(c_ice, snow_lambda(tt),...
                precip_out(:,tt), snow_mwe(tt), snow_e(tt), air_T_out(tt), ro_water,...
                T_melt, ro_snow_initial, snow_ro(tt - 1), snow_mwe(tt - 1), k_ice,...
                ro_ice, b_exp, c_water, snow_l_ro(tt), k_snow_ice_max);
        end


        %surface flux based on snow occurance
        if track_b(2) == 0
            [q(tt)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                T(1,tt), ice_albedo, AWS_albedo, Io_ice, e_ice,0);
        elseif track_b(2) == 1
            [snow_albedo] = snow_albedo_calc(snow_albedo, precip_out(tt),...
                refresh_albedo, t_step, new_snow_albedo, snow_albedo_min, tau_snow_cold,...
                tau_snow_melt, snow_T(tt), T_melt);
            snow_albedo = snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end

            [q(tt)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                snow_T(tt), snow_albedo, AWS_albedo, Io_ice, e_ice,0);
        end

        %apply surface flux to surface cell
        if track_b(2) == 0
            enthalpy(1,tt) = enthalpy(1,tt) + q(tt)*t_step*3600; %*t_step... to get energy over timestep
        elseif track_b(2) == 1
            snow_e(tt) = snow_e(tt) + q(tt)*t_step*3600; %*t_step... to get energy over timestep
        else; fprintf('Error in main, stage 4, surface flux application. \n \n');
        end


        %calcualte temperature and lambda for cells based on enthalpy
        [T(:,tt), lambda(:,tt), model_stage, snow_T(tt)] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water, snow_mwe(tt), snow_e(tt),lake_ind);

        if model_stage == 5; track_a(tt) = 5; track_b(2) = 0; end %to leapfrog temp_lambda_profile that shoots it back to model_stage = 4

        %calculate lake index for conduction
        [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
            slush_lid_bot] = lake_index_new(lambda(:,tt), ice_lim, total_grid_num, grid_profile(:,tt),...
            tt, t_step, slush_lid_threshold, slush_lid_num);

        %record depth and average temperature
        lake_depth(tt) = bot_lake_depth;
        lake_T_av(tt) = sum(T(lake_ind,tt).*grid_profile(lake_ind,tt))./sum(grid_profile(lake_ind,tt));
        lid_thick(tt) = sum(grid_profile(1:min(lake_ind) - 1,tt));

        %calculate turbulent heat flux if required
        if lake_depth(tt) >= lake_turb_threshold && lake_depth(tt) < lake_prof_threshold
            lake_mode = 1; %mode 1 = turbulent flux only, no convection

            [enthalpy(:,tt)] = turbulent_flux1(lake_ind, lake_T_av(tt), air_T_out(tt),...
                t_step, ro_water, c_water, J, T_melt, lake_depth(tt), Lf,...
                enthalpy(:,tt), grid_profile(:,tt), T(:,tt), c_ice, lake_mode, model_stage);


            %recalcualte temperature and lambda for cells based on enthalpy
            [T(:,tt), lambda(:,tt), model_stage, snow_T(tt)] = temp_lambda_profile1(enthalpy(:,tt),...
                total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
                ro_water, c_water, snow_mwe(tt), snow_e(tt),lake_ind);

        elseif lake_depth(tt) >= lake_prof_threshold %calculate convection profile if required
            lake_mode = 2; %mode 2 = turbulent flux and convection

            %calculate turbulent heat flux
            [~, dT, lf_lower(tt), lf_upper(tt)] = turbulent_flux1(lake_ind, lake_T_av(tt), air_T_out(tt),...
                t_step, ro_water, c_water, J, T_melt, lake_depth(tt), Lf,...
                enthalpy(:,tt), grid_profile(:,tt), T(:,tt), c_ice, lake_mode, model_stage);

            %apply temperature change to lake
            T(lake_ind,tt) = T(lake_ind,tt) + dT;

            %calculate convection profile and assign to main temperature and
            %enthalpy profiles
            [T(lake_ind,tt), enthalpy(lake_ind,tt)] = convection1(T(lake_ind,tt),...
                grid_profile(lake_ind,tt), ies80, ro_water, c_ice, T_melt, Lf, c_water);

            %apply upper and lower convective flux
            enthalpy(max(lake_ind) + 1,tt) = enthalpy(max(lake_ind) + 1,tt) + lf_lower(tt);
            if min(lake_ind) - 1 >= 1 %to prevent if no lid occurs just before stage 4 switch
                enthalpy(min(lake_ind) - 1,tt) = enthalpy(min(lake_ind) - 1,tt) + lf_upper(tt);
            end

            %recalcualte temperature and lambda for cells based on enthalpy
            [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
                total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
                ro_water, c_water, snow_mwe(tt), snow_e(tt),lake_ind);

        end

        if track_a(tt) == 5; model_stage = 5;  track_b(2) = 0; snow_z(tt) = 0; end %complete leapfrog from initial temp_lambda_profile

        %set thermal conductivity and heat capacity of ice cells
        k(:,tt) = lambda(:,tt)*k_water + (1 - lambda(:,tt))*k_ice;
        c(:,tt) = lambda(:,tt)*c_water + (1 - lambda(:,tt))*c_ice;

        %conduction calculations. Model stage held at 4 to prevent a
        %repetitive initialisation section at the start of stage 5
        [enthalpy_temp] = conduct_update_enthalpy_new(k(:,tt), c(:,tt), T(:,tt),...
            grid_profile(:,tt), enthalpy(:,tt), t_step, total_grid_num, ro_water, lower_boundary,...
            4,  snow_T(tt), snow_k(tt), snow_c(tt), snow_z(tt), snow_e(tt));

        %seperate enthalpy_temp and apply snow updates if required
        if track_b(2) == 0
            enthalpy(:,tt + 1) = enthalpy_temp;
        elseif track_b(2) == 1 %so snow is occuring
            snow_e(tt) = enthalpy_temp(1);
            enthalpy(:,tt + 1) = enthalpy_temp(2:end);

            %second round of snow updates
            [snow_ro(tt), snow_z(tt), snow_lambda(tt), snow_T(tt), snow_ice(tt), snow_l_ro(tt)] =...
                update_snow_b(snow_e(tt), T_melt, snow_mwe(tt), c_ice, ro_water,...
                Lf, snow_ro(tt), ro_melt_max, t_step, tau_ro, ro_cold_max,...
                snow_ice(tt-1), snow_lambda(tt-1), snow_l_ro(tt));

            %snow updates
            snow_T(tt + 1) = snow_T(tt); %for surface flux calculations
            snow_e(:,tt + 1) = snow_e(:,tt);
            snow_lambda(:,tt + 1) = snow_lambda(:,tt);
            snow_mwe(:,tt + 1) = snow_mwe(:,tt);
            snow_ro(:,tt + 1) = snow_ro(:,tt);
            snow_l_ro(:,tt + 1) = snow_l_ro(:,tt);
            snow_z(:,tt + 1) = snow_z(:,tt);
            snow_ice(:,tt + 1) = snow_ice(:,tt);

        else; fprintf('Error in main, stage 4 enthalpy allocation. \n \n');
        end

        %updates
        grid_profile(:,tt + 1) = grid_profile(:,tt);
        T(1,tt + 1) = T(1,tt); %only top cell gets updated for surface flux calculations
        track_a(1,tt + 1) = 4; %model stage = 2
        

        if model_stage == 5
            fprintf('Moving from stage 4 to %1.0f. Day %3.0f. \n \n', model_stage, tt*t_step/24)

            %reset snow trackers
            track_b(1:2) = 0;
        end
    end

    %STAGE 5, lid breakup
    if model_stage == 5

        %calculate tau profile
        [tau(:,tt)] = tau_calc(tau_water, tau_slush, tau_ice, total_grid_num, lambda(:,tt));

        %calculate lake albedo
        [lake_albedo(tt)] = cal_lakealbedo(albedo_num,lake_depth(tt));
        %albedo multiplier
        lake_albedo(tt) = lake_albedo(tt)*albedo_mult; if lake_albedo(tt) > 1; lake_albedo(tt) = 1; end
        
        if gcnet_num == 1 %then calculate LW_in seperately
            [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1,tt), model_stage);
        end
        %surface flux
        [q(tt)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
            wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
            T(1,tt), lake_albedo(tt), 0, Io_water, e_ice,3); %AWS_albedo = 0 as lake used

        %calculate Io_profile
        [Io(:,tt)] = Io_calc(Io_water, Io_ice, Io_slush, total_grid_num, lambda(:,tt));

        %SW propagation through water. 0 is to prevent AWS albedo from
        %being used
        [SW_prop(:,tt)] = SW_propagate1(SW_down_out(tt), grid_profile(:,tt), albedo_out(tt),...
            t_step, Io(:,tt), tau(:,tt), lake_albedo(tt), 0, basal_SW_distribute);

        %apply incoming energy flux to enthalpy array
        enthalpy(1,tt) = enthalpy(1,tt) + q(tt)*3600*t_step;
        enthalpy(:,tt) = enthalpy(:,tt) + SW_prop(:,tt);

        %calcualte temperature and lambda for cells based on enthalpy
         [T(:,tt), lambda(:,tt), model_stage, snow_T(tt)] = temp_lambda_profile1(enthalpy(:,tt),...
            total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
            ro_water, c_water, snow_mwe(tt), snow_e(tt),lake_ind);

        %calculate lake index for conduction
        [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
            slush_lid_bot] = lake_index_new(lambda(:,tt), ice_lim, total_grid_num, grid_profile(:,tt),...
            tt, t_step, slush_lid_threshold, slush_lid_num);

        %see if lid breakup should occur if it hasn't already
        if track_b(3) == 0
            [track_b(3)] = lid_instability(lambda(:,tt), grid_profile(:,tt), slush_lid_top,...
                slush_lid_bot, slush_mech_threshold, breakup_threshold, ro_water, ro_ice);
        end

        %lid may also breakup if lake_index subfunction combines lakes into
        %one. If this is the case, prompt lid breakup and reobtain lake
        %index with a greater slush_lid_threshold
        if (min(lake_ind) - 1) <= 1
            track_b(3) = 1;

            [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                slush_lid_bot] = lake_index_new(lambda(:,tt), ice_lim, total_grid_num, grid_profile(:,tt),...
                tt, t_step, 0.5, slush_lid_num);
        end

        %record lake values
        surface_lake_depth(tt) = top_lake_depth; %depth of the surface lake
        lake_depth(tt) = bot_lake_depth; %depth of the main lake
        sur_lake_bot = max(top_lake_ind) + 1; %bottom of the surface lake
        lake_top = min(lake_ind) - 1; %top of the main lake
        lake_bot = max(lake_ind) + 1; %bottom of the main lake

        %if lid is unstable
        if track_b(3) == 1

            %lid breakup and lake combine
            [model_stage, lake_ind, enthalpy(:,tt + 1), track_b, track_a(1,tt + 1),...
                lake_depth(tt)] = lid_breakup(track_b, enthalpy(:,tt),...
                lake_ind, grid_profile(:,tt));

            %updates
            fprintf('Moving from stage 5 to 3. Day %3.0f. \n \n', tt*t_step/24)
            T(:,tt + 1) = T(:,tt); %skip temperature calculation and just update

        else %calculate average temperatures as normal

            sur_lake_T_av(tt) = sum(grid_profile(top_lake_ind,tt).*T(top_lake_ind,tt))./...
                sum(grid_profile(top_lake_ind,tt)); %surface lake average temperature
            lake_T_av(tt) = sum(grid_profile(lake_ind,tt).*T(lake_ind,tt))./...
                sum(grid_profile(lake_ind,tt)); %main lake average temperature
        end

        snow_z(tt) = 0;
        snow_mwe(tt) = 0;
        %only continue if lid has not yet disintegrated
        if model_stage == 5 || model_stage == 4

            %calculate turbulent heat flux for surface lake if required.
            %convection profile is not considered here
            if surface_lake_depth(tt) >= lake_turb_threshold

                %calculate upper and lower turbulent heat fluxes following
                %Buzzard (2017)
                if min(top_lake_ind) == 1
                    lf_upper(tt) = sign(sur_lake_T_av(tt) - air_T_out(tt))*ro_water*c_water*J*abs(sur_lake_T_av(tt) - air_T_out(tt))^(4/3)*3600*t_step;
                else
                    lf_upper(tt) = sign(sur_lake_T_av(tt) - T_melt)*ro_water*c_water*J*abs(sur_lake_T_av(tt) - T_melt)^(4/3)*3600*t_step;
                end
                lf_lower(tt) = sign(sur_lake_T_av(tt) - T(max(lake_ind)+1,tt))*ro_water*c_water*J*abs(sur_lake_T_av(tt) - T(max(lake_ind)+1,tt))^(4/3)*3600*t_step; %assumes ice-water interface at T_melt

                %calculate change in temperature of core
                dT = (-lf_upper(tt) - lf_lower(tt))/...
                    (ro_water*c_water*surface_lake_depth(tt));

                %apply turbulent fluxes as required
                if min(top_lake_ind) > 1 %if statement to apply top flux if lid has formed
                    enthalpy(min(top_lake_ind) - 1,tt) = enthalpy(min(top_lake_ind) - 1,tt) + lf_upper(tt);
                end
                enthalpy(max(top_lake_ind)+1,tt) = enthalpy(max(top_lake_ind)+1,tt) + lf_lower(tt);

                %calcualte temperature and lambda for cells based on enthalpy
                [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
                    total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe(tt), snow_e(tt),lake_ind);

                %update core temperature
                sur_lake_T_av(tt) = sur_lake_T_av(tt) + dT;
                T(top_lake_ind,tt) = sur_lake_T_av(tt);

            end

            %calculate turbulent heat flux for main lake if required
            if lake_depth(tt) >= lake_turb_threshold && lake_depth(tt) < lake_prof_threshold
                lake_mode = 1; %mode 1 = turbulent flux only, no convection

                [enthalpy(:,tt)] = turbulent_flux1(lake_ind, lake_T_av(tt), air_T_out(tt),...
                    t_step, ro_water, c_water, J, T_melt, lake_depth(tt), Lf,...
                    enthalpy(:,tt), grid_profile(:,tt), T(:,tt), c_ice, lake_mode, model_stage);

                %recalcualte temperature and lambda for cells based on enthalpy
                [T(:,tt), lambda(:,tt), model_stage] = temp_lambda_profile1(enthalpy(:,tt),...
                    total_grid_num, grid_profile(:,tt), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe(tt),snow_e(tt),lake_ind);


            elseif lake_depth(tt) >= lake_prof_threshold %calculate convection profile for main lake if required
                lake_mode = 2; %mode 2 = turbulent flux and convection

                %calculate turbulent heat flux
                [~, dT, lf_lower(tt)] = turbulent_flux1(lake_ind, lake_T_av(tt), air_T_out(tt),...
                    t_step, ro_water, c_water, J, T_melt, lake_depth(tt), Lf,...
                    enthalpy(:,tt), grid_profile(:,tt), T(:,tt), c_ice, lake_mode, model_stage);

                %apply temperature change to lake
                T(lake_ind,tt) = T(lake_ind,tt) + dT;


                %calculate convection profile and assign to main temperature and
                %enthalpy profiles
                [T(lake_ind,tt), enthalpy(lake_ind,tt)] = convection1(T(lake_ind,tt),...
                grid_profile(lake_ind,tt), ies80, ro_water, c_ice, T_melt, Lf, c_water);

                %apply lower convective flux
                enthalpy(max(lake_ind) + 1,tt) = enthalpy(max(lake_ind) + 1,tt) + lf_lower(tt);

                %a recalc_ind is used as if updating the whole profile the
                %old enthalpy updates the new surface temperature
                recalc_ind = [min(lake_ind) - 1; lake_ind; max(lake_ind + 1)];
                if model_stage == 4; track_a(tt) = 4; end %leapfrog stage 4 to avoid reset in temp_lambda_profile due to recalc_ind

                %recalcualte temperature and lambda for lake cells based on enthalpy
                [T(recalc_ind,tt), lambda(recalc_ind,tt), model_stage] = temp_lambda_profile1(enthalpy(recalc_ind,tt),...
                    numel(recalc_ind), grid_profile(recalc_ind,tt), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe(tt), snow_e(tt),lake_ind);

                if track_a(tt) == 4; model_stage = 4; end %copmlete leapfrog

            end

            %set thermal conductivity and heat capacity of ice cells
            k(:,tt) = lambda(:,tt)*k_water + (1 - lambda(:,tt))*k_ice;
            c(:,tt) = lambda(:,tt)*c_water + (1 - lambda(:,tt))*c_ice;
            if model_stage == 4
                track_b(2) = 0;
            end
            %conduction calculations
            snowz = 0;
            [enthalpy(:,tt+1)] = conduct_update_enthalpy_new(k(:,tt), c(:,tt), T(:,tt),...
                grid_profile(:,tt), enthalpy(:,tt), t_step, total_grid_num, ro_water, lower_boundary,...
                model_stage, snow_T(tt), snow_k(tt), snow_c(tt), snowz, snow_e(tt));
        end

        %updates
        grid_profile(:,tt + 1) = grid_profile(:,tt);
        T(1,tt + 1) = T(1,tt); %only top cell gets updated for surface flux calculations
        lambda(:,tt + 1) = lambda(:,tt);
        track_a(1,tt + 1) = model_stage;

    end
end

toc; fprintf('\n \n')

%% 
bm = sum(lambda(19:end,:).*grid_profile(19:end,:));%cumulative bottom ablation
basal_freeze = max(bm)-bm(t_num-1);

% %FIGURES 
% 
% %set empty values as nans
% empty_ind = T == 0;
% enthalpy(empty_ind) = nan;
% lambda(empty_ind) = nan;
% SW_prop(empty_ind) = nan;

%conversions. This gets very messy due to Matlab not having functionality
%for multiple colourmaps on one figure. Best not to fiddle with this (except
%for the bits which are explicitly fiddlable)!
%
T_c = T - T_melt; %convert from K to C for readability 
T_c(T_c < 0 & T_c >= -0.004) = -0.004; %basically, it seems the resolution in zero threshold clips temperatures below ~0.004 and above 0 as water, not slush, so just correcting this for the plot
T_c(1,1) = 11;
T_c(T_c < low_plot_T) = low_plot_T; %remove values lower than threshold, allows splitting into ice and water but prevents observing erros as easily
T_c(T_c > hi_plot_T) = hi_plot_T;
T_c_ice = T_c; T_c_lake = T_c; T_c_temp = T_c; %ice, lake, and temporary arrays
T_c_temp(T == 0) = max(max(T_c)) + 2; %so that area before hydrograph is light grey
% 
T_c_temp(lambda < 0.3 & lambda > 0) = max(max(T_c)) - (1/5000)*abs(diff([low_plot_T max(max(T_c))]));
T_c_temp(lambda < 1 & lambda > 0) = max(max(T_c)) - (1/5000)*abs(diff([low_plot_T max(max(T_c))])); %to force the value of 0 into the second lowest colour of the colour map
T_c_range = diff([max(max(T_c)) min(min(T_c))]); %temperature range, range doesn't work with 2D arrays
zero_threshold = abs(low_plot_T/(T_c_range)); %fraction value of where 0 is
T_c_ice(lambda > 0) = nan; %omit slush and water cells from plot
T_c_lake(lambda == 0) = nan; %omit ice and slush cells from plot
T_c_lake(T == 0) = max(max(T_c)) + 2; %so that area before hydrograph is light grey
% 
T_c_lake(lambda < 0.3 & lambda > 0) = max(max(T_c)) - (1/5000)*abs(diff([low_plot_T max(max(T_c))])); %to force the value of 0 into the second lowest colour of the colour map
T_c_lake(lambda < 1 & lambda > 0) = max(max(T_c)) - (1/5000)*abs(diff([low_plot_T max(max(T_c))])); %to force the value of 0 into the second lowest colour of the colour map
T_c_slush = lambda < 0.3 & lambda > 0; %value of one where slush exists 
% 
T_c_slush = lambda < 1 & lambda > 0; %value of one where slush exists 
lake_z_plot = lake_depth + sub_lake_depth;
lake_z_plot(lake_z_plot == 0) = nan;
lid_thick(lid_thick == 0) = nan;
lid_thick = lid_thick - min(lid_thick); %bring to surface
y_labx = -25; %percent for plotting, needs fiddling


if output_figs == 1

    %create colourmap for main figure. This only works if lake has formed, i.e.
    %if T_c is above zero at any point
    if sum(any(T_c > 0)) > 0
        jet_sec = jet(12000);
        jet_sec = jet_sec(1:10000,:);
        gray_sec = gray(5000);
        gray_sec = gray_sec(2000:4300,:);
        ice_map = jet_sec(2:floor(zero_threshold*10000),:); %the 2 is to allow black as the first colour of lake_map
        water_map = gray_sec(1:ceil((1 - zero_threshold)*10000) - 2,:); %the -1 is to allow light grey as last colour
        lake_map = [0.15 0.15 0.15; ice_map; water_map; 0.33 0.33 0.33; 1 1 1];
        custom_colormap = 1; %to prevent another search of T_c
    else
        custom_colormap = 0;
    end

%     create axis
    x_axis = (0:t_step:t_total)/24; %time in days
    x_axis = x_axis(1:end - 1);
    y_axis = ice_grid_z:ice_grid_z:((total_grid_num - deep_ice_grid_num)*ice_grid_z);%linspace(grid_profile(1,t_num - 1), sum(grid_profile(1:ice_grid_num - deep_ice_grid_num,t_num - 1)), (ice_grid_num - deep_ice_grid_num)); %depth in metres
    if rem(t_total, t_step) ~= 0
         disp('t_step error, does not fit into t_total')
    end

%     main figure
    figure (5)
    fs = 8; %fontsize
set (gcf,'Unit','centimeters','position',[8,8,15,12])
    %flot surface energy flux
    sp1 = subplot(4, 1, 1);
    sp1a = plot(x_axis, q);
    set(sp1a, 'color', [0.8 0.8 0.8]);
    set(sp1a, 'linewidth', 0.05);
    set(sp1, 'Position', [0.1,0.79,0.8,0.20]);
    set(sp1, 'color', [1 1 1])
    %set(gca, 'XTick',[]);
    set(gca,'XTickLabel',[])
    xlim([min(x_axis) max(x_axis)])
    ylim([(min(q) - abs(min(q)*0.2)) (max(q) + max(q)*0.2)])
    hold on
    mov_q = movmean(q, 24/t_step);
    sp1b = plot(x_axis, mov_q, 'k');
    set(sp1b, 'linewidth', 0.8);
    sp1_lab = ylabel(sp1, 'Energy Flux (W/m^2)');
    set(sp1_lab, 'FontSize', fs)
%     pos = get(sp1_lab, 'Pos');
%     set(sp1_lab, 'Pos', [y_labx pos(2) pos(3)])
    set(sp1_lab, 'Pos',[-39.34860191317142,-99.0230650053864,-1])
    set(gca,'FontSize',fs)

    if sum(any(snow_z > 0)) > 0 %only plot if snow has occured
        sp2 = subplot(4, 1, 2);
        snow_z(snow_z == 0) = nan;
        sp2a = plot(x_axis, snow_z, 'k'); %total depth of snow
        xlim([min(x_axis) max(x_axis)])
        ylim([0 1.05*max(snow_z)])
        sp2_lab = ylabel(sp2, 'Snow Depth (m)');
        set(sp2_lab, 'FontSize', fs)
%         pos = get(sp2_lab, 'Pos'); % INTI CODE
%         set(sp2_lab, 'Pos', [y_labx pos(2) pos(3)])
        set(gca,'yaxislocation','right');
%         set(sp2_lab, 'Pos',[-32.17430095658568,0.27505174539488,-1])
        hold on
        snow_water_plot = snow_mwe.*snow_lambda;
        snow_water_ind = snow_water_plot > 0; %index for addition
        snow_water_plot(snow_water_ind) = snow_water_plot(snow_water_ind);
        snow_water_plot(snow_water_plot == 0) = nan;
        sp2c = plot(x_axis, snow_water_plot, '--b');
        sp2_labx = xlabel(sp2, 'Time (days)');
        set(sp2_labx, 'FontSize', fs)
        xlim([min(x_axis) max(x_axis)])
        hold off
        set(gca,'XTickLabel',[])
%         set(gca, 'XTick',[]);
        %set(sp2, 'Position', [0.1 0.68 0.8 0.10]);
        set(sp2,'Position',[0.1,0.61,0.8,0.17])
        set(sp2, 'color', [1 1 1])
        set(gca,'FontSize',fs)
    end

    sp3 = subplot(4, 1, 3);
    if custom_colormap == 1
        colormap(lake_map)
    else
        colormap jet 
    end
    img1 = imagesc(x_axis, y_axis, T_c_temp(1:(total_grid_num - deep_ice_grid_num),:));
    xlim([min(x_axis) max(x_axis)])
    caxis([low_plot_T max(max(T_c))])
    sp3_labx = xlabel(sp3, 'Time (days)');
    set(sp3_labx, 'FontSize', fs)
    sp3_laby = ylabel(sp3, 'Depth (m)');
    set(sp3_laby, 'FontSize', fs)
    pos = get(sp3_laby, 'Pos');
    set(sp3_laby, 'Pos', [y_labx pos(2) pos(3)])
%     set(sp3, 'Position', [0.1 0.03 0.8 0.59]);
    set(sp3, 'Position',[0.1,0.186,0.8,0.413])
    T_bar = colorbar('southoutside');
    s = sprintf('Temperature (%cC)', char(176));
    set(get(T_bar,'XLabel'), 'String', s)
    set(get(T_bar,'XLabel'), 'FontSize', fs)
    set(T_bar, 'Position',[0.1,0.08,0.8,0.03])
    set(gca,'FontSize',fs)
    hold on 
end
%%
figure(7)
basal_ablation = bm;
%load ablation rate
ablation_10_ted = readtable('ablation_10.txt');
ablation_10_ted = table2array(ablation_10_ted);
% ablation_10_ted = ablation_10_ted(3:end,:);
con1 = ablation_10_ted(1,2);
ablation_10_ted(:,2) = ablation_10_ted(:,2)-con1;

diffs = abs(x_axis - ablation_10_ted(1,1));

[~, idx] = min(diffs);
con = ablation_10_ted(1,2)-basal_ablation(idx);
basal_ablation = basal_ablation+con;
hold on
plot(ablation_10_ted(:,1),ablation_10_ted(:,2),'.-')
plot(x_axis,basal_ablation,'.-')
xlim([175,180])
xlabel('day')
ylabel('Bottom ablation (m)')

bluesnow_m = interp1(x_axis,basal_ablation,ablation_10_ted(:,1));
prec = ablation_10_ted(:,2);
R2m = 1-(sum((bluesnow_m-prec).^2)./sum((prec-mean(prec)).^2));
rmsem = sqrt(mean((bluesnow_m-prec).^2));
biasm = mean(bluesnow_m-prec);

% if albedo_num == 0
%     save('BSr2a0.mat')
% else
%     save('BSr2a2.mat')
% end
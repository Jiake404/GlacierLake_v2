function [runoff,melt,q,qm,SW_flux] = runoff_prod(SW_down_out, LW_down_out, air_T_out, rh_out, hum_out, pressure_out, wind_speed_out,...
    albedo_out, precip_out, run_day, t_step, albedo_mult, model_stage,  ini_condtions, ice_grid_num,...
    ice_grid_z, deep_ice_grid_num, deep_ice_grid_z, lower_boundary, snow_threshold, new_snow_albedo, snow_albedo_min,...
    refresh_albedo, tau_snow_cold, tau_snow_melt, AWS_albedo, ice_albedo, Io_ice, Io_water, Io_slush, ro_snow_initial,...
    b_exp, ro_melt_max, ro_cold_max, tau_ro, k_snow_ice_max, T_melt, Lf, ro_water, c_ice, c_water, k_water, k_ice,...
    k_air, e_ice, ro_ice, gcnet_num, AWS)
tic

% 
t_total = run_day*24;   %days to hours
t_num = t_total/t_step; %total number of time steps
total_grid_num = ice_grid_num + deep_ice_grid_num;  %total number of grid cells
grid_profile = zeros(total_grid_num, t_num);        %(m) array of depths of each cell in grid
grid_profile(1:ice_grid_num,1) = ice_grid_z;
grid_profile(ice_grid_num + 1:total_grid_num,1) = deep_ice_grid_z;

%INITIALISE
%create arrays
enthalpy = zeros(total_grid_num, t_num);        %(J) enthalpy of cells
lambda = zeros(total_grid_num, t_num);          %water content (0 = no water, 1 = all water)
T = zeros(total_grid_num, t_num);               %(K) temperature
k = zeros(total_grid_num, t_num);               %(W/(m.K)) thermal conductivity
c = zeros(total_grid_num, t_num);               %(J/K) heat capacity
q = zeros(1, t_num);                            %(J) energy transfer at the surface
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
track_a = zeros(1, t_num);                      %tracking array, keep everything in this one and add rows as required.

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

%row 1 = model stage
%initialise arrays
T(:,1) = linspace(ice_T_surface, ice_T_bottom, total_grid_num)';
enthalpy(:,1) = c_ice.*ro_water.*T(:,1).*grid_profile(:,1);    %following Benedek (2014), all ice to begin with as default
lambda(:,1) = 0;                                                    %begin as ice as default
track_a(1,1) = model_stage;                                           %to avoid 0 error


% create arrays
qm = zeros(1,t_num);
latm = zeros(1,t_num);
melt = zeros(1,t_num);
con = zeros(1,t_num);
liquid1 = zeros(1,t_num);
q_max = zeros(1,t_num);
grid_num = zeros(1,t_num);
runoff = zeros(1,t_num);
snow_mass = zeros(1,t_num);
w_mi = zeros(1,t_num);
w_mwe = zeros(1,t_num);
dw = zeros(1,t_num);
dh = zeros(1,t_num);
ice_in_s = zeros(1,t_num);
mwe_add = zeros(1,t_num);
%warnings
if Io_ice > Io_water || Io_slush > Io_water; fprintf('Warning, edit SW_prop to allow for a lower value of Io_water than Io_ice. \n \n'); end

%MAIN LOOP
fprintf('Main loop commencing. \n \n')
%% 
for tt = 1:t_num-1

    pre = precip_out(:,tt);
    SWD = SW_down_out(tt);
    AL = albedo_out(tt);
    LWD = LW_down_out(tt);
    ta = air_T_out(tt);
    RH = rh_out(tt);
    HUM = hum_out(tt);
    P = pressure_out(tt);
    WS = wind_speed_out(tt);

    snow_mwe(tt) = snow_mwe(tt)+pre(2,1);
    % switch to model stage 2 if threshold has been passed
    if snow_mwe(tt) >= snow_threshold
        model_stage = 2;
    else
        model_stage = 1;
    end

    %STAGE 1, bare ice
    if model_stage == 1
         if gcnet_num == 1 %then calculate LW_in seperately
            [LWD] = longwave_in1(RH, ta);
         end
        %surface flux
         [q(tt),sens_flux(tt,1),lat_flux(tt,1),SW_flux(tt,1),LW_down(tt,1),LW_out(tt,1),CT(tt,1)]...
             = surface_flux1(SWD, AL, ta, WS, P, LWD, HUM, T(grid_num(tt)+1,tt), ice_albedo, 1, Io_ice, e_ice,0);

        if q(tt) > 0
           [qm(tt),sens_fluxm(tt,1),lat_fluxm(tt,1),SW_fluxm(tt,1),LW_downm(tt,1),LW_outm(tt,1)] = surface_flux1(SWD,...
               AL, ta, WS, P,LWD, HUM, T_melt, ice_albedo, 1, Io_ice, e_ice,0);
        end

        if qm(tt) > 0
            melt(tt) = qm(tt)*3600/(Lf*ro_water);
            if latm(tt)<0
                con(tt) = abs(latm(tt))*3600/(2.5*10^9);
            end
        end

        liquid1(tt) = pre(3)+con(tt);
        % The maximum of heat released when liquid water is completely frozen
        q_max(tt) = liquid1(tt)*ro_water*Lf/(3600);

        % Assume all liquid water freezes
        % apply incoming energy flux to enthalpy array
        enthalpy(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt) + (q(tt)+q_max(tt))*3600*t_step;
        
        %calcualte temperature for cells based on enthalpy
        if enthalpy(grid_num(tt)+1,tt)<0 || grid_profile(grid_num(tt)+1,tt)<0.01
           
            enthalpy(grid_num(tt)+2,tt) = enthalpy(grid_num(tt)+2,tt)+enthalpy(grid_num(tt)+1,tt);
            
            grid_profile(grid_num(tt)+2,tt) = grid_profile(grid_num(tt)+2,tt)+grid_profile(grid_num(tt)+1,tt)+liquid1(tt)*ro_water/ro_ice;% 无径流
            
            grid_num(tt) = grid_num(tt)+1;
            
            T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
        else 
           
            grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)+liquid1(tt)*ro_water/ro_ice;% 无径流
            
            T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
        end

        % The first grid does not reach the melting point, no runoff
        if T(grid_num(tt)+1,tt) <= T_melt
            runoff(tt) = 0;
        end
        
        % The first grid reaches the melting point and produces runoff
        if T(grid_num(tt)+1,tt) > T_melt
            
            runoff(tt) = (enthalpy(grid_num(tt)+1,tt)-T_melt*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt))/(Lf*ro_water);% m w.e.
            
            dh(tt) = runoff(tt)*ro_water/ro_ice; 
            
            if dh(tt) < grid_profile(grid_num(tt)+1,tt)
                
                T(grid_num(tt)+1,tt) = T_melt;
                
                grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)-dh(tt);
                
                enthalpy(grid_num(tt)+1,tt) = T(grid_num(tt)+1,tt)*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt);
            else 
                
                enthalpy(grid_num(tt)+2,tt) = enthalpy(grid_num(tt)+2,tt)+enthalpy(grid_num(tt)+1,tt);
                
                grid_profile(grid_num(tt)+2,tt) = grid_profile(grid_num(tt)+2,tt)+grid_profile(grid_num(tt)+1,tt);
                
                grid_num(tt) = grid_num(tt)+1;
                
                T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)./(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
                
                if T(grid_num(tt)+1,tt) > T_melt
                    
                    T(grid_num(tt)+1,tt) = T_melt;
                    
                    runoff(tt) = (enthalpy(grid_num(tt)+1,tt)-T_melt*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt))/(Lf*ro_water);% m w.e.
                    
                    dh(tt) = runoff(tt)*ro_water/ro_ice; % m
                    
                    grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)-dh(tt);
                    
                    enthalpy(grid_num(tt)+1,tt) = T_melt*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt);
                    
                    if grid_profile(grid_num(tt)+1,tt)<0
                        fprintf('spin_up may require 3 cells to merge, which is not considered')
                        disp(tt)
                        break
                    end
                else
                    runoff(tt) = 0;
                end
            end
        end
        
        if grid_num(tt)>0
            enthalpy(1:grid_num(tt),tt) = nan;
            T(1:grid_num(tt),tt) = nan;
            grid_profile(1:grid_num(tt),tt) = 0;
        end

        %update thermal conductivity and heat capacity to that of ice
        lambda(:,tt) = 0;
        k(:,tt) = k_ice;
        c(:,tt) = c_ice;

        %conduction calculations
        [enthalpy(grid_num(tt)+1:total_grid_num,tt+1)] = conduct_update_enthalpy1(k(grid_num(tt)+1:total_grid_num,tt), c(grid_num(tt)+1:total_grid_num,tt),...
            T(grid_num(tt)+1:total_grid_num,tt),grid_profile(grid_num(tt)+1:total_grid_num,tt), enthalpy(grid_num(tt)+1:total_grid_num,tt), ...
            t_step, total_grid_num-grid_num(tt), ro_water, lower_boundary, model_stage);
       
        if grid_num(tt)>0
            enthalpy(1:grid_num(tt),tt+1) = nan;
            k(1:grid_num(tt),tt+1) = nan;
            c(1:grid_num(tt),tt+1) = nan;
            T(1:grid_num(tt),tt+1) = nan;
        end
        
        %updates
        grid_profile(:,tt + 1) = grid_profile(:,tt);
        T(:,tt + 1) = T(:,tt); %only top cell gets updated for surface flux calculations
        snow_mwe(1,tt + 1) = snow_mwe(:,tt);
        grid_num(tt+1) = grid_num(:,tt);
        track_a(1,tt + 1) = 1; %model stage = 1
    end

    %STAGE 2, snow on ice
    if model_stage == 2
        
        if track_a(1,tt) ~= 2 %if coming from another stage then initialise

            %initialise snow variables
            [snow_lambda(tt), snow_ro(tt), snow_c(tt), snow_T(tt),...
                snow_e(tt), snow_k(tt), snow_z(tt), snow_ice(tt), snow_l_ro(tt)] =...
                initialise_snow1(snow_mwe(tt),ro_snow_initial,...
                c_ice, k_ice, k_air, ro_water, T_melt, ta);
            snow_albedo = new_snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end

        elseif track_a(1,tt) == 2 %if within stage 2 then update as normal
            if tt == 1
                [snow_c(tt), snow_e(tt), snow_lambda(tt), snow_ro(tt),...
                    snow_k(tt)] = update_snow_a1(c_ice, snow_lambda(tt),...
                    pre, snow_mwe(tt), snow_e(tt), ta, ro_water,...
                    T_melt, ro_snow_initial, snow_ro(tt), snow_mwe(tt), k_ice,...
                    ro_ice, b_exp, c_water, snow_l_ro(tt), k_snow_ice_max);
            else
                %update snow variables
                [snow_c(tt), snow_e(tt), snow_lambda(tt), snow_ro(tt),...
                    snow_k(tt)] = update_snow_a1(c_ice, snow_lambda(tt),...
                    pre, snow_mwe(tt), snow_e(tt), ta, ro_water,...
                    T_melt, ro_snow_initial, snow_ro(tt - 1), snow_mwe(tt - 1), k_ice,...
                    ro_ice, b_exp, c_water, snow_l_ro(tt), k_snow_ice_max);
            end

        end

        %surface flux 
        % snow albedo
        [snow_albedo] = snow_albedo_calc(snow_albedo, pre(2),...
            refresh_albedo, t_step, new_snow_albedo, snow_albedo_min, tau_snow_cold,...
            tau_snow_melt, snow_T(tt), T_melt);
        snow_albedo = snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end

        % calculate LW_in
        if gcnet_num == 1 %then calculate LW_in seperately
            [LWD] = longwave_in1(RH, ta);
        end


        [q(tt),sens_flux(tt,1) , lat_flux(tt,1) , SW_flux(tt,1) , LW_down(tt,1) ,...
            LW_out(tt,1),CT(tt,1)] = surface_flux1(SWD, AL, ta, WS, P,...
            LWD, HUM, snow_T(tt), snow_albedo, 1, Io_ice, e_ice,0);

        if q(tt) > 0
            [qm(tt),sens_fluxm(tt,1) , lat_fluxm(tt,1) , SW_fluxm(tt,1) , LW_downm(tt,1) ,...
                LW_outm(tt,1)] = surface_flux1(SWD, AL, ta, WS, P,...
                LWD, HUM, T_melt, ice_albedo, 1, Io_ice, e_ice,0);
        end

        if qm(tt) > 0
            melt(tt) = qm(tt)*3600/(Lf*ro_water); %m w.e.

            if latm(tt)<0
                con(tt) = abs(latm(tt))*3600/(2.5*10^9);
            end
        end

        liquid1(tt) = pre(3)+con(tt);
        q_max(tt) = liquid1(tt)*ro_water*Lf;

        % Assume all liquid water freezes
        snow_e(tt) = snow_e(tt)+q_max(tt)+q(tt)*3600;
        snow_T(tt) = snow_e(tt)/(c_ice*ro_water*snow_mwe(tt));

        if grid_profile(grid_num(tt)+1,tt)<0.01
            enthalpy(grid_num(tt)+2,tt) = enthalpy(grid_num(tt)+2,tt)+enthalpy(grid_num(tt)+1,tt);
            grid_profile(grid_num(tt)+2,tt) = grid_profile(grid_num(tt)+2,tt)+grid_profile(grid_num(tt)+1,tt)+liquid1(tt)*ro_water/ro_ice;% 无径流
            grid_num(tt) = grid_num(tt)+1;
            T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
        end

        if snow_T(tt) <= 273.15&&liquid1(tt)>0
            runoff(tt) = 0;
            % The maximum height of ice that can be frozen in a snow layer
            h_ice = snow_z(tt)*((ro_ice-snow_l_ro(tt))/ro_ice);

            % The height at which liquid water is completely frozen
            h_water = ro_water*liquid1(tt)/ro_ice;
            
            if h_water >= h_ice
                enthalpy(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)+c_ice*ro_water*snow_T(tt)*(snow_z(tt)+h_water-h_ice);
                grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)+snow_z(tt)+h_water-h_ice;
                T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
                snow_z(tt) = 0;
                snow_mwe(tt) = 0;
                snow_e(tt) = 0;
                model_stage = 1;
            else
                % The actual height of liquid water frozen in snow
                ice_in_s(tt) = h_water*ro_ice/ro_water;

            end
        end

        if snow_T(tt) > 273.15
            runoff(tt) = (snow_e(tt)-c_ice*ro_water*T_melt*snow_mwe(tt))./(ro_water*Lf);%mwe
            if runoff(tt) <= liquid1(tt) 
                % If runoff < liquid1, indicating that a portion of liquid1 freezes and releases heat 
                % so that the first cell reaches the melting point and then produces runoff
                h_water = (liquid1(tt)-runoff(tt))*ro_water/ro_ice;

                % The maximum height of ice that can be frozen in a snow layer
                h_ice = snow_z(tt)*((ro_ice-snow_l_ro(tt))/ro_ice);
                if h_water >= h_ice
                    enthalpy(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)+c_ice*ro_water*T_melt*(snow_z(tt)+h_water-h_ice);
                    grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)+snow_z(tt)+h_water-h_ice;
                    T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
                    snow_z(tt) = 0;
                    snow_mwe(tt) = 0;
                    snow_e(tt) = 0;
                    model_stage = 1;
                else
                    ice_in_s(tt) = h_water*ro_ice/ro_water;                    
                end

            elseif runoff(tt)>liquid1(tt)&&runoff(tt)<liquid1(tt)+snow_mwe(tt)
                snow_T(tt) = T_melt;
                snow_mwe(tt) = snow_mwe(tt)-(runoff(tt)-liquid1(tt));
                if snow_mwe(tt)<snow_threshold
                    snow_z(tt) = 0;
                    snow_mwe(tt) = 0;
                    snow_e(tt) = 0;
                    model_stage = 1;
                else
                    snow_e(tt) = snow_T(tt)*c_ice*ro_water*snow_mwe(tt);
                    snow_z(tt) = snow_mwe(tt)*(1 - snow_lambda(tt))*(ro_water/snow_l_ro(tt));
                end

            elseif runoff(tt) >= liquid1(tt)+snow_mwe(tt)
                % If this condition is met, the snow layer is completely melted
                runoff_res = runoff(tt)-liquid1(tt)-snow_mwe(tt);

                % The remaining energy after melting the snow layer.
                en_res = runoff_res*ro_water*Lf;

                % whether the remaining energy can bring the first ice grid to the melting point.
                enthalpy(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)+en_res;
                
                T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
                
                if T(grid_num(tt)+1,tt) <= T_melt
                    runoff(tt) = liquid1(tt)+snow_mwe(tt);
               
                else
                   
                    drf = (enthalpy(grid_num(tt)+1,tt)-T_melt*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt))/(Lf*ro_water);
                    
                    dh(tt) = drf*ro_water/ro_ice;
                    
                    if dh(tt) < grid_profile(grid_num(tt)+1,tt)
                        
                        runoff(tt) = liquid1(tt)+snow_mwe(tt)+drf;% m w.e.
                        
                        grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)-dh(tt);
                        
                        T(grid_num(tt)+1,tt) = T_melt;
                        
                        enthalpy(grid_num(tt)+1,tt) = T_melt*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt);
                    else
                        
                        runoff(tt) = liquid1(tt)+snow_mwe(tt)+grid_profile(grid_num(tt)+1,tt);% m w.e.
                        
                        dh1 = (dh(tt)-grid_profile(grid_num(tt)+1,tt))*ro_ice/ro_water;
                        
                        de = dh1*Lf*ro_water;
                        
                        grid_num(tt) = grid_num(tt)+1;
                        
                        enthalpy(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)+de;
                        
                        T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
                        
                        if T(grid_num(tt)+1,tt) > T_melt
                            
                            runoff(tt) = runoff(tt)+(enthalpy(grid_num(tt)+1,tt)-T_melt*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt))/(Lf*ro_water);% m w.e.
                            
                            dh(tt) = ((enthalpy(grid_num(tt)+1,tt)-T_melt*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt))/(Lf*ro_water))*ro_water/ro_ice; %由于产流损失的格子高度
                            
                            if dh(tt) < grid_profile(grid_num(tt)+1,tt)
                                
                                T(grid_num(tt)+1,tt) = T_melt;
                                
                                grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)-dh(tt);
                                
                                enthalpy(grid_num(tt)+1,tt) = T(grid_num(tt)+1,tt)*c_ice*ro_water*grid_profile(grid_num(tt)+1,tt);
                            else 
                                fprintf('spin_up may require 3 cells to merge, which is not considered')
                                disp(tt)
                            end
                        end
                    end
                end
                snow_z(tt) = 0;
                snow_mwe(tt) = 0;
                snow_e(tt) = 0;
                model_stage = 1;
            end
        end
       
        % Calculate meltwater retention
        if model_stage == 2 && runoff(tt) > 0
            % Calculate snow layer mass
            snow_mass(tt) = snow_z(tt)*snow_l_ro(tt);%kg

            %Calculate the maximun retention in snow layer from Coléou and Lesaffre (1998)
            w_mi(tt) = 0.017+0.057*(ro_ice/snow_l_ro(tt)-1);
            w_mwe(tt) = snow_mass(tt).*w_mi(tt)/ro_water; %retent water, m w.e.

            % Calculate the actual runoff and the water retention in the snow layer
            if track_a(1,tt) ~= 2
                if runoff(tt) > w_mwe(tt)
                    runoff(tt) = runoff(tt)-w_mwe(tt);
                    dw(tt) = w_mwe(tt);% actual water retention
                else
                    runoff(tt) = 0;
                    dw(tt) = runoff(tt); 
                end
                snow_mwe(tt) = snow_mwe(tt)+dw(tt);

            elseif track_a(1,tt) == 2 

                if w_mwe(tt) > dw(tt) && runoff(tt) >= w_mwe(tt)-dw(tt)
                    runoff(tt) = runoff(tt)-(w_mwe(tt)-dw(tt));
                    mwe_add(tt) = w_mwe(tt)-dw(tt);
                    dw(tt) = w_mwe(tt);
                    
                elseif w_mwe(tt) > dw(tt) && runoff(tt) < w_mwe(tt)-dw(tt)
                    mwe_add(tt) = runoff(tt);
                    dw(tt) = dw(tt)+runoff(tt);
                    runoff(tt) = 0;
                    
                elseif w_mwe(tt) <= dw(tt)
                    runoff(tt) = runoff(tt)+(dw(tt)-w_mwe(tt));
                    mwe_add(tt) = -(dw(tt)-w_mwe(tt));
                    dw(tt) = w_mwe(tt);
                    
                end
            end
        end

        snow_mwe(tt) = snow_mwe(tt)+ice_in_s(tt)+mwe_add(tt);
        snow_e(tt) = snow_T(tt)*c_ice*ro_water*snow_mwe(tt);
        snow_T(tt) = snow_e(tt)/(c_ice*ro_water*snow_mwe(tt));
        snow_z(tt) = snow_mwe(tt)*(1 - snow_lambda(tt))*(ro_water/snow_l_ro(tt));
        snow_ice(tt) = ice_in_s(tt)/snow_mwe(tt);
        snow_lambda(tt) = dw(tt)/snow_mwe(tt);
        
        if snow_ice(tt)>=1
            model_stage = 1;
            enthalpy(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)+c_ice*ro_water*snow_T(tt)*snow_z(tt);
            grid_profile(grid_num(tt)+1,tt) = grid_profile(grid_num(tt)+1,tt)+snow_z(tt);
            T(grid_num(tt)+1,tt) = enthalpy(grid_num(tt)+1,tt)/(c_ice*ro_water*grid_profile(grid_num(tt)+1,tt));
        end
        
        lambda(:,tt) = 0;
        T(grid_num(tt)+1:total_grid_num,tt) = enthalpy(grid_num(tt)+1:total_grid_num,tt)...
            ./(c_ice.*ro_water.*grid_profile(grid_num(tt)+1:total_grid_num,tt));

        %set thermal conductivity and heat capacity of ice cells
        k(grid_num(tt)+1:total_grid_num,tt) = lambda(grid_num(tt)+1:total_grid_num,tt)*...
            k_water + (1 - lambda(grid_num(tt)+1:total_grid_num,tt))*k_ice;
        c(grid_num(tt)+1:total_grid_num,tt) = lambda(grid_num(tt)+1:total_grid_num,tt)*c_water...
            + (1 - lambda(grid_num(tt)+1:total_grid_num,tt))*c_ice;

        %seperate enthalpy_temp
        if model_stage == 2
%             tt
            %conduction calculations
            [enthalpy_temp] = conduct_update_enthalpy1(k(grid_num(tt)+1:total_grid_num,tt), c(grid_num(tt)+1:total_grid_num,tt), T(grid_num(tt)+1:total_grid_num,tt),...
                grid_profile(grid_num(tt)+1:total_grid_num,tt), enthalpy(grid_num(tt)+1:total_grid_num,tt), t_step, total_grid_num-grid_num(tt), ro_water, lower_boundary,...
                model_stage, snow_T(tt), snow_k(tt), snow_c(tt), snow_z(tt), snow_e(tt));
            snow_e(tt + 1) = enthalpy_temp(1);
            enthalpy(grid_num(tt)+1:total_grid_num,tt + 1) = enthalpy_temp(2:total_grid_num-grid_num(tt)+1);

            [snow_ro(tt), snow_z(tt), snow_T(tt), snow_l_ro(tt)] =...
                update_snow_b2(snow_lambda(tt), snow_ice(tt), snow_e(tt), T_melt, snow_mwe(tt), c_ice, ro_water,...
                snow_ro(tt), ro_melt_max, t_step, tau_ro, ro_cold_max, snow_l_ro(tt));

            if grid_num(tt)>0
                enthalpy(1:grid_num(tt),tt) = nan;
                T(1:grid_num(tt),tt) = nan;
                grid_profile(1:grid_num(tt),tt) = 0;
            end

            %updates
            dw(:,tt+1) = dw(tt);
            grid_profile(:,tt + 1) = grid_profile(:,tt);
            grid_num(:,tt+1) = grid_num(:,tt);
            T(:,tt + 1) = T(:,tt); %only top cell gets updated for surface flux calculations
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
            track_a(1,tt + 1) = 1; %model stage = 1
            %conduction calculations
            [enthalpy(grid_num(tt)+1:total_grid_num,tt+1)] = conduct_update_enthalpy1(k(grid_num(tt)+1:total_grid_num,tt), c(grid_num(tt)+1:total_grid_num,tt),...
                T(grid_num(tt)+1:total_grid_num,tt),grid_profile(grid_num(tt)+1:total_grid_num,tt), enthalpy(grid_num(tt)+1:total_grid_num,tt), ...
                t_step, total_grid_num-grid_num(tt), ro_water, lower_boundary, model_stage);
            if grid_num(tt)>0
                enthalpy(1:grid_num(tt),tt+1) = nan;
                k(1:grid_num(tt),tt+1) = nan;
                c(1:grid_num(tt),tt+1) = nan;
                T(1:grid_num(tt),tt+1) = nan;
            end
            %updates
            grid_profile(:,tt + 1) = grid_profile(:,tt);
            T(:,tt + 1) = T(:,tt); %only top cell gets updated for surface flux calculations
            snow_mwe(1,tt + 1) = snow_mwe(:,tt);
            grid_num(tt+1) = grid_num(:,tt);
            track_a(1,tt + 1) = 1; %model stage = 1
        end
    end
end
%%
% sum(melt)
% sum(runoff)
% toc; fprintf('\n \n')
% for i = 1:size(runoff,2)/24
%     melt_d(i,1) = sum(melt(i*24-23:i*24));
%     runoff_d(i,1) = sum(runoff(i*24-23:i*24));
% end
% filename = [AWS,'runoffd3.mat'];
% save(filename,'runoff_d');
end
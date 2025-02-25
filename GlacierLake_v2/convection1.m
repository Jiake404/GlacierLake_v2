function [lake_T_prof, enthalpy_out] = convection1(T, grid_profile, ies80, ro_water, c_ice, T_melt, Lf, c_water)
%GlacierLake function to apply UNESCO IES80 algorithm to lake temperature
%profile following Saloranta and Andersen (2007) 
    
    %create lake temperature index  (convec_T) 

    convec_T = T;
    convec_z = grid_profile;

    %calculate density profile in kg/m3
    convec_ro = polyval(ies80, max(0, convec_T(:))) + min(convec_T(:), 0);

    %calculate how much lighter a layer is than the one below. My lake treats
    %the last cell as positive to account for dense sediment. Field
    %observations of melt pools suggests although slush is less dense, it
    %has some structural integrity that keeps it from immediately moving to
    %the surface.
    d_ro = [diff(convec_ro); 1];

    %while loop to reshuffle layers
    while any(d_ro < 0)

        unstable_layers = d_ro <= 0; %1 = layer is heaver or equal to layer below
        unstable_start = find(diff([0; unstable_layers]) == 1); %index of where unstable column(s) start
        unstable_end = find(diff([0; unstable_layers]) == -1) - 1; %index of where unstable column(s) end

        %average the density of unstable sections of the temperature profile
        for i = 1:length(unstable_start)
            unstable_ind = [unstable_start(i):unstable_end(i) + 1]; %index of cells to be used
            T_mix = sum(convec_T(unstable_ind).*convec_z(unstable_ind))./...
                sum(convec_z(unstable_ind)); %average the temperature of these cells
            convec_T(unstable_ind) = ones(numel(unstable_ind),1)*T_mix; %apply this temperature back to the cells
        end

        %recalculate density
        convec_ro = polyval(ies80, max(0, convec_T(:))) + min(convec_T(:), 0);
        d_ro = [diff(convec_ro); 1];
    end

    %spring/autumn overturn is omitted as all water would need to rise to
    %3.9 deg C
    
    %output
    lake_T_prof = convec_T;
    
    %recalculate enthalpy
    enthalpy_out = ro_water.*convec_z.*...
    (c_ice*T_melt + Lf + c_water.*lake_T_prof - c_water*T_melt);










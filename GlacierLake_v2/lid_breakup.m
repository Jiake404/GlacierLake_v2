function [model_stage, lake_ind, enthalpy, track_b, track_a, lake_depth] =...
    lid_breakup(track_b, enthalpy, lake_ind, grid_profile)
%GlacierLake function to break up lid and combine surface and main lakes

    %calculate the new lake index which is now a surface lake
    lake_ind = 1:max(lake_ind);

    %calculate the new average temperature. This should at least be
    %above freezing unless there are any ice cells left
    lake_enthalpy_total = sum(enthalpy(lake_ind));

    %calculate the new lake depth
    lake_depth = sum(grid_profile(lake_ind)); %+ ...
        %grid_profile(sur_lake_bot,tt)*lambda(sur_lake_bot,tt);

    %assign new lake enthalpy according to cell thickness
    enthalpy(lake_ind) = lake_enthalpy_total*...
        (grid_profile(lake_ind)/sum(grid_profile(lake_ind)));

    %go back into stage 3
    model_stage = 3;
    track_b(3) = 0;
    track_a = 5; %set stage tracker as 5



















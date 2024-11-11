function [lid_breakup] = lid_instability(lambda, grid_profile, slush_lid_top,...
    slush_lid_bot, slush_mech_threshold, breakup_threshold, ro_water, ro_ice)
%GlacierLake function to determine stability/instability of lid

    %obtain logical arrays of lid strength based on value of
    %slush_mech_threshold
    slush_lid_strong_ind = lambda <= slush_mech_threshold;

    %if there are 'strong' slush cells within the lid
    if sum(slush_lid_strong_ind(1:slush_lid_bot)) >= 1

        %find base of strong section
        slush_lid_strong_bot = slush_lid_bot - find(flipud(slush_lid_strong_ind(1:slush_lid_bot)) == 1, 1, 'first') + 1;

        %find top of strong section
        slush_lid_strong_top = slush_lid_strong_bot - find(flipud(slush_lid_strong_ind(1:slush_lid_strong_bot)) == 0, 1, 'first') + 2;

        %so base of weak section is one above top of strong section
        slush_lid_weak_bot = slush_lid_strong_top - 1;

        %calculate density of weak section weighted by depth
        ro_weak = sum((lambda(1:slush_lid_weak_bot).*ro_water + (1 - lambda(1:slush_lid_weak_bot)).*ro_ice).*...
            grid_profile(1:slush_lid_weak_bot))./sum(grid_profile(1:slush_lid_weak_bot));

        %calculate density of strong section of lid
        ro_strong = sum((lambda(slush_lid_strong_top:slush_lid_strong_bot).*ro_water + (1 - lambda(slush_lid_strong_top:slush_lid_strong_bot)).*ro_ice).*...
            grid_profile(slush_lid_strong_top:slush_lid_strong_bot))./sum(grid_profile(slush_lid_strong_top:slush_lid_strong_bot));

        %is the weak section above the strong section heavier
        weight_weak = ro_weak*sum(grid_profile(1:slush_lid_weak_bot));
        weight_strong = ro_strong*sum(grid_profile(slush_lid_strong_top:slush_lid_strong_bot));
        
        %and then a relationship between the two can be established
        if weight_strong/weight_weak >= breakup_threshold

            %set tracker
            lid_breakup = 1;

        else

            %set tracker
            lid_breakup = 0;
        end
        
    %if slush_lid_top is emtpy then there is no slush lid to breakup    
    elseif isempty(slush_lid_top)
        
        %set tracker
        lid_breakup = 0;
        
    %if there are no strong cells within the lid then it should
    %disintegrate
    elseif ~ismember(1, slush_lid_strong_ind(1:slush_lid_bot))

        %set tracker
        lid_breakup = 1;
    end
end
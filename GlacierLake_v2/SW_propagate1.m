function [SW_prop] = SW_propagate1(SW_down, grid_profile,...
    albedo_out, t_step, Io, tau, albedo, AWS_albedo, basal_SW_distribute)
%GlacierLake subfunction to calculate SW propagation through ice or
%water-ice column

    %set albedo to AWS albedo if required
    if AWS_albedo == 1
        albedo = albedo_out;
    end
    
    cum_depth = cumsum(grid_profile);
    pc = find(abs(diff(Io))>0); %do any phase changes exist, if so where
    SW_prop = zeros(numel(grid_profile),1);

    if isempty(pc) %if no phase changes (i.e. no water)
        SW_in = Io(1)*(1- albedo)*SW_down*t_step.*3600; %SW_in per time step
        SW_prop_total = SW_in.*exp(tau.*-cum_depth); %Beer Lambert law
    else %if water is present
        
        %initialisations
        SW_prop_vector = zeros(numel(grid_profile),numel(pc) + 1);
        SW_in = zeros(1,numel(pc) + 1);
        SW_boundary_ab = zeros(1,numel(pc)); %this only applys to the phase change within the column whereas SW_in also applys to the surface
        SW_prop_total = zeros(numel(grid_profile),1);
        SW_in(1) = Io(1)*(1- albedo)*SW_down*t_step.*3600;% 此处为什么是用albedo_out而不是albedo
        SW_prop_vector(:,1) = SW_in(1).*exp(tau.*-cum_depth); 

        %loop to account for phase changes
        for i = 1:numel(pc)

            %find SW propagating just before phase change, using the
            %entire amount to make it so far
            SW_in(i + 1) = Io(pc(i) + 1)*SW_prop_vector(pc(i),i); % - SW_prop_vector(pc(i) + 1,i));
            SW_boundary_ab(i) = (1 - Io(pc(i) + 1))*SW_prop_vector(pc(i),i); % - SW_prop_vector(pc(i) + 1,i));

            %calculate SW propagation for next phase
            SW_prop_vector(pc(i) + 1:end,i + 1) = SW_in(i + 1).*exp(tau(pc(i) + 1:end).*-cum_depth(pc(i) + 1:end));
        end

        SW_prop_total(1:pc(1)) = SW_prop_vector(1:pc(1),1);

        %update pc to include end if necessary for next step. This
        %forces the following for loop to fill the enitre array.
        pc_all = pc;
        if pc(end) ~= numel(grid_profile)
            pc_all(end + 1) = numel(grid_profile);
        else
            fprintf('If this is happening the boundary condition in SW_propagate subfunction has gone really really wrong, investigate. \n \n')
        end

        %combine SW_prop_vector
        for i = 1:numel(pc_all) - 1
            SW_prop_total(pc_all(i) + 1:pc_all(i + 1)) = ...
                SW_prop_vector(pc_all(i) + 1:pc_all(i + 1),i + 1);
        end
        
    end
    
    %subtract SW to obtain SW at each grid cell. End value is ignored which
    %is fine as it's a constant boundary condition. If basal_SW_distribute
    %= 1 then SW goes into water just above lake bottom, rather than the
    %ice/slush below
    if basal_SW_distribute == 1
        %this means that the SW will go into water just above lake
        %bottom rather than ice/slush just below
        SW_prop(1:end - 1) = SW_prop_total(1:end - 1) - SW_prop_total(2:end);
    else
        SW_prop(2:end) = SW_prop_total(1:end - 1) - SW_prop_total(2:end);
    end
    
end
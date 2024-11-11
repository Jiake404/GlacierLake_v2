function [tau] = tau_calc(tau_water, tau_slush, tau_ice, total_grid_num, lambda)
%GlacierLake function to calculate tau profile

    tau = zeros(total_grid_num,1);
    tau(lambda >= 1) = tau_water;
    tau(lambda < 1 & lambda > 0) = tau_slush;
    tau(lambda == 0) = tau_ice;
    
end

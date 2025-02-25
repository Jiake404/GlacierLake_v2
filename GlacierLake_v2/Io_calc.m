function [Io] = Io_calc(Io_water, Io_ice, Io_slush, total_grid_num, lambda)
%GlacierLake function to calculate Io profile

    Io = zeros(total_grid_num,1);
    Io(lambda >= 1) = Io_water;
    Io(lambda < 1 & lambda > 0) = Io_slush;
    Io(lambda == 0) = Io_ice;
    
end

%function to caluclate snow albedo, following Essery (2015)

function [snow_albedo] = snow_albedo_calc(snow_albedo, precip_out,...
    refresh_albedo, t_step, snow_albedo_max, snow_albedo_min, tau_snow_cold,...
    tau_snow_melt, surface_T, T_melt)

%calculate albedo decay timescale
if surface_T < T_melt
    tau = tau_snow_cold*3600; %*3600 to get seconds
else
    tau = tau_snow_melt*3600;
end

%calculate rt, reciprocal timescale for albedo adjustment
rt = 1/tau + (precip_out/3600)/refresh_albedo; %/3600 to get seconds

%calculate limiting albedo
alim = (snow_albedo_min/tau + ((precip_out/3600)*snow_albedo_max)/refresh_albedo)/rt;

%calculate albedo
snow_albedo = alim + (snow_albedo - alim)*exp(-rt*t_step*3600); %*3600 to get seconds

%apply adjustments if range has been exceeded
if snow_albedo < snow_albedo_min; snow_albedo = snow_albedo_min; end
if snow_albedo > snow_albedo_max; snow_albedo = snow_albedo_max; end

end
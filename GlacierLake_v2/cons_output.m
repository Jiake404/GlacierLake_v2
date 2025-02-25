function [t_step, model_stage, albedo_mult, precip_mult, albedo_add, ice_grid_num, ice_grid_z,...
    deep_ice_grid_num, deep_ice_grid_z, lower_boundary, low_plot_T, hi_plot_T, ice_lim,...
    snow_threshold, slush_lid_num, new_snow_albedo, snow_albedo_min, max_hydro_input,...
    breakup_threshold, lake_turb_threshold, lake_prof_threshold, slush_lid_threshold,...
    slush_mech_threshold, stage_print_hysteresis, refresh_albedo, tau_snow_cold, tau_snow_melt,...
    output_figs, sp_year, AWS_albedo, basal_SW_distribute,...
    ice_albedo, Io_ice, Io_water, Io_slush, tau_ice, tau_water, tau_slush, tau_snow, J, ro_snow_initial,...
    b_exp, ro_melt_max, ro_cold_max, tau_ro, k_snow_ice_max, T_melt, Lf, ro_water, c_ice, c_water,...
    k_water, k_ice, k_air, e_ice, T_ro_max, ro_ice, ies80, t_total, t_num, total_grid_num, grid_profile,...
    lake_grid_num,hydro_T,ice_T_bottom,ice_T_surface] = cons_output(run_day)

t_step = 1;           %(hrs) time step
ice_T_bottom = -5;      %(deg C) basal ice temperature 
ice_T_surface = -4;     %(deg C) initial ice temperature 
model_stage = 1;        %model initially set at 1 as default
albedo_mult = 1;        %multiple for increasing or decreasing albedo. Currently no check on wether this goes above 1
precip_mult = 1;        %precipitation multiplier 
albedo_add = 0;         %addition or subtraction for AWS albedo data
ice_grid_num = 150;     %nuber of ice cells
ice_grid_z = 0.1;       %(m) thickness of ice cells
deep_ice_grid_num = 30; %number of deep ice cells at coarser resolution
deep_ice_grid_z = 1;    %(m) thickness of deep ice cells
lower_boundary = 265;   %(K) lower boundary condition
low_plot_T = -35;       %(deg C) lowest temp to display in figures
hi_plot_T = 10;         %(deg C) highest temp to display in figures
ice_lim = 1;            %wubbold, not totally sure what this does, test at end
snow_threshold = 0.01;  %(mwe snow) depth at which snow is incorporated 
slush_lid_num = 3;      %
new_snow_albedo = 0.75; %albedo of fresh snow
snow_albedo_min = 0.5;  %minimum snow albedo
max_hydro_input = 5;    %(m) maximum depth of water to be put into the hydrograph function
breakup_threshold = 0.2;% how much heaver area above mechanically strong lid can be before failure
lake_turb_threshold = 0.09;     %(m) depth above which turbulent heat flux is considered
lake_prof_threshold = 1.6;     %(m) depth where profiling is coducted
slush_lid_threshold = 0.01;     %(m) if lid is entirely made of slush, this is the cut off for what is considered to be part of the lid
slush_mech_threshold = 0.01;    %proportion of water permisable in slush cell before it is considered mechanically weak
stage_print_hysteresis = 20;    %(days) how many days must pass before stage switch is displayed (in case it switches back and forth multiple times) 
refresh_albedo = 10;            %(kg/m2) snowfall required to refresh albedo
tau_snow_cold = 1000;           %(hr) cold snow albedo decay timescale
tau_snow_melt = 100;            %(hr) melting snow albedo decay timescale
sp_year = 5;                    %(yr) the number of spin-up years

%USER SWITCHES
output_figs = 1;        %if = 1 then show model output figures, if = 0 do not
AWS_albedo = 0;         %if = 1 then use albedo from AWS, if = 0 use from user input albedo
basal_SW_distribute = 1;%if = 1 then distribute the SW that reaches the ice-water interface in stage 3     

%PARAMETERS
ice_albedo = 0.55;
Io_ice = 0;             %fraction of SW not absorbed at the surface for Beer Lambert law
Io_water = 0.6;        %fraction of SW not absorbed at the surface for Beer Lambert law
Io_slush = Io_water;    %fraction of SW not absorbed at the surface for Beer Lambert law
tau_ice = 1.5;          %(1/m) bulk shortwave extinction coefficient
tau_water = 0.025;      %(1/m) bulk shortwave extinction coefficient
tau_slush = tau_water;  %(1/m) bulk shortwave extinction coefficient
tau_snow = 1.4;         %(1/m) bulk shortwave extinction coefficient. 1.4 from Grenfell and Maykut (1977). Varies between 1.1 - 1.5
J = 1.907*10^-5;        %(m/sK^(1/3)) turbulent heat flux factor
ro_snow_initial = 180;  %(kg/m^3) this can vary a lot, from Essery FSM
b_exp = 2;              %thermal conductivity exponent,from Essery FSM
ro_melt_max = 500;      %(kg/m^3) maximum density if snow is at T_melt
ro_cold_max = 300;      %(kg/m^3) maximum density if snow is below T_melt
tau_ro = 20*3600;       %(s) compaction timescale of snow
hydro_T = 0.8;          %(deg C) temperature of incoming water in hydrograph
k_snow_ice_max = 1.5;   %(W/mK) maximum allowable thermal conductivity of snow ice

%CONSTANTS
T_melt = 273.15;        %(K) melt temperature of water
Lf = 3.348*10^5;        %(J/kg) latent heat of fusion of water
ro_water = 1000;        %(kg/m^3) density of water 
c_ice = 2108;           %(J/kgK) specific heat capacity of ice 
c_water = 4217;         %(J/kgK) Specific heat capacity of water 
k_water = 0.569;        %(W/mK) conductivity of water 
k_ice = 1.88;           %(W/mK) conductivity of ice
k_air = 0.022;          %(W/mK) conductivity of air, from Buzzard (2017) following Moaveni (2010)
e_ice = 0.99;           %emissivity of bare ice
T_ro_max = 277.13;      %(K) temperature of maximum water density
ro_ice = 900;           %(kg/m^3) density of ice

%International Equation of State (UNESCO, 1981). 5-order polynomial for 
%density as function of temperature. Following MyLake
ies80 = [6.536332e-9, -1.120083e-6, 1.001685e-4, -9.09529e-3, 6.793952e-2, 999.842594];

%INITIAL CONVERSIONS, CALCULATIONS, AND INITIALISATIONS
t_total = run_day.*24;   %days to hours
t_num = t_total/t_step; %total number of time steps
total_grid_num = ice_grid_num + deep_ice_grid_num;  %total number of grid cells
ice_T_surface = ice_T_surface + T_melt;             %convert to kelvin
ice_T_bottom = ice_T_bottom + T_melt;               %convert to kelvin
grid_profile = zeros(total_grid_num, t_num);        %(m) array of depths of each cell in grid
grid_profile(1:ice_grid_num,1) = ice_grid_z; 
grid_profile(ice_grid_num + 1:total_grid_num,1) = deep_ice_grid_z; 
lake_grid_num = max_hydro_input/ice_grid_z;         %the max number of cells that the hydrograph will use
% resize_tt = floor((resize_day*24/t_step));          %time step for resizing array to incorporate hydrograph input
hydro_T = hydro_T + T_melt;                         %convert to kelvin
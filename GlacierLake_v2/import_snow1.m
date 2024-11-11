function [precip_out,precip_s] = import_snow1(AWS,clip_day1,clip_day2,run_start,run_end,t_num)
%GlacierLake subfunction to import and plot snow and precipitation data from RACMO NetCDF files. Data
%in comes in daily values. If this is not the case edits must be made.
str1 = 'SwissCamp';
if strcmp(AWS, str1)
    AWS = 'Bluesnow';
end

run_starty = year(datetime(run_start));
run_endy = year(datetime(run_end));
yn = run_starty:run_endy;

snowfall_data = [];
precip_data = [];
for i = 1:length(yn)
    % Precipitation and snowfall data should be prepared in advance
    eval(strcat("load('precip_",AWS,"_",num2str(yn(i)),".mat');")) % mm w.e. per day
    eval(strcat("load('snowfall_",AWS,"_",num2str(yn(i)),".mat');")) % mm w.e. per day
    eval(strcat("precip_data = [precip_data;precip_",AWS,"_",num2str(yn(i)),"];"))
    eval(strcat("snowfall_data = [snowfall_data;snowfall_",AWS,"_",num2str(yn(i)),"];"))
end
%mm w.e. to m w.e.
precip_data = precip_data./1000; % m w.e.per day
snowfall_data = snowfall_data./1000; % m w.e.per day

%get rid of negative values
precip_data(precip_data<0) = 0;
snowfall_data(snowfall_data<0) = 0;

% 
prec_out = (precip_data(clip_day1:clip_day2))';
snow_out = (snowfall_data(clip_day1:clip_day2))';
precip_out = [prec_out;snow_out;prec_out-snow_out];

%resize to fit time step used in main.m
original_length = length(snow_out);
precip_out = imresize(precip_out, [3 t_num], 'bilinear'); %reshape to the correct size
precip_out = precip_out*(original_length/t_num); %correct for stretching
precip_out(precip_out < 0) = 0; %remove negative precipitation values

%spin-up data
clip = yeardays(run_start(1))*24;
precip_s = precip_out(1:3,1:clip);










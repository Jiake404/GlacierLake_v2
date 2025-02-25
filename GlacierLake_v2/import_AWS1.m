function [SW_down_out, LW_down_out, air_T_out, rh_out, hum_out, pressure_out,...
    wind_speed_out, albedo_out, s_data, time_out] = import_AWS1(AWS,clip_hour1,clip_hour2,run_start,gcnet_num)

%GlacierLake subfunction to import data from PROMICE weather data and plot if required 

%turn off warning
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

%load data here and add extra if required. PROMICE data available from: 
%https://www.promice.dk/WeatherStations.html
Data = readtable([AWS '_hour.csv']); %1 
time = year(Data{:,1});
% Data = table2array(Data(:,2:end));
% A=zeros(size(Data,1),4);
% Data=[A,Data];

if gcnet_num == 0
    %air temperature
    air_T = Data{:,'t_u'};
    air_T(air_T == 0) = NaN;
    air_T = air_T + 273.15; %convert to kelvin
    %remove -999s through linear interpolation
    temp_var = air_T;
    ind = 1:length(temp_var);
    ix = logical(temp_var > -473.15); %-200 degrees
    temp_var2 = temp_var;
    air_T=interp1(ind(ix),temp_var2(ix),ind,'linear');
    air_T = air_T';

    %SW down
    SW_down = Data{:,'dsr'};
    SW_down(SW_down < 0) = 0;
    %remove -999s through linear interpolation
    temp_var = SW_down;
    ind = 1:length(temp_var);
    ix = logical(temp_var >= 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    SW_down=interp1(ind(ix),temp_var2(ix),ind,'linear');
    SW_down=SW_down';

    %SW out
    SW_out = Data{:,'usr'};
    SW_out(SW_out < 0) = 0;
    %remove -999s through linear interpolation
    temp_var = SW_out;
    ind = 1:length(temp_var);
    ix = logical(temp_var >= 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    SW_out=interp1(ind(ix),temp_var2(ix),ind,'linear');
    SW_out=SW_out';  

    % albedo
    albedo = SW_out./SW_down;
%     albedo = Data{:,'albedo'};
    %remove -999s through linear interpolation
    temp_var = albedo;
    ind = 1:length(temp_var);
    ix = logical(temp_var > 0 & temp_var < 1); %0 to 1 range
    temp_var2 = temp_var;
    albedo=interp1(ind(ix),temp_var2(ix),ind,'linear');
    albedo=albedo';

    %LW down
    LW_down = Data{:,'dlr'};
    LW_down(LW_down == 0) = NaN;
    %remove -999s through linear interpolation
    temp_var = LW_down;
    ind = 1:length(temp_var);
    ix = logical(temp_var > 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    LW_down=interp1(ind(ix),temp_var2(ix),ind,'linear');
    LW_down=LW_down';

%     %LW out
%     LW_out = Data{:,'ulr'};
%     LW_out(LW_out == 0) = NaN;
%     %remove -999s through linear interpolation
%     temp_var = LW_out;
%     ind = 1:length(temp_var);
%     ix = logical(temp_var > 0); %-200 degrees chosen but may need to change, data missing = -999
%     temp_var2 = temp_var;
%     LW_out=interp1(ind(ix),temp_var2(ix),ind,'linear');
%     LW_out=LW_out';
%     LW_out(LW_out>316)=316;
%     Ts = (LW_out./(5.67*10.^(-8))).^0.25;


    %pressure
    pressure = Data{:,'p_u'};
    pressure(pressure == 0) = NaN;
    pressure = pressure/10; %convert from hPa to kPa
    %remove -999s through linear interpolation
    temp_var = pressure;
    ind = 1:length(temp_var);
    ix = logical(temp_var > 0); %0 chosen as physically impossible
    temp_var2 = temp_var;
    pressure=interp1(ind(ix),temp_var2(ix),ind,'linear');
    pressure=pressure';

    %relative humidity
    rh = Data{:,'rh_u_cor'};
    rh(rh == 0) = NaN;
    %remove -999s through linear interpolation
    temp_var = rh;
    ind = 1:length(temp_var);
    ix = logical(temp_var > -473.15); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    rh=interp1(ind(ix),temp_var2(ix),ind,'linear');
    rh=rh';

    %convert relative humidity to specific humidity following https://earthscience.stackexchange.com/questions/5076/how-to-calculate-specific-humidity-with-relative-humidity-temperature-and-pres
    %constants
    Rv = 461.5; %J/kg*K specific gas constant
    Lv = 2257*10^3; %J/kg specific enthalpy of vaporisation
    T0 = 273.15; %K reference temperature, may need changing
    Rd = 287.058; %J/kg*K specific gas constant for dry air
    es0 = 613; %Pa water saturation pressure at 0 degrees

    %rh humidity from % to proportion
    rh_proportion = rh./100;

    %loop
    es = es0.*exp((Lv./Rv)*(1./T0 - 1./air_T));
    e = rh_proportion.*es;
    w = (e.*Rd)./(Rv.*(pressure.*1000 - e)); %*1000 for Kpa to pa
    q = w./(w + 1);

    q(q == 0) = NaN; %get rid of zeros that have crept in
    hum = q; %save specific humidity as hum in kg/kg

    %windspeed
    wind_speed = Data{:,'wspd_u'};
    wind_speed(wind_speed == 0) = NaN;
    %remove -999s through linear interpolation
    temp_var = wind_speed;
    ind = 1:length(temp_var);
    ix = logical(temp_var > 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    wind_speed=interp1(ind(ix),temp_var2(ix),ind,'linear');
    wind_speed=wind_speed';

    %
    SW_down_out = SW_down(clip_hour1:clip_hour2);
    SW_down_out(SW_down_out < 0) = 0;
    albedo_out = albedo(clip_hour1:clip_hour2);
    albedo_out(isnan(albedo_out)) = 0;
    LW_down_out = LW_down(clip_hour1:clip_hour2);
    air_T_out = air_T(clip_hour1:clip_hour2);
    rh_out = rh(clip_hour1:clip_hour2);
    hum_out = hum(clip_hour1:clip_hour2);
    pressure_out = pressure(clip_hour1:clip_hour2);
    wind_speed_out = wind_speed(clip_hour1:clip_hour2);
%     LW_up_out = LW_out(clip_hour1:clip_hour2);
%     Ts_out = Ts(clip_hour1:clip_hour2);
%     Ts_out(Ts_out>273.15)=273.15;
    time_out = time(clip_hour1:clip_hour2);
else
    %air temperature
    air_T = Data{:,'TA1'};
    air_T(air_T == 0) = NaN;
    air_T = air_T + 273.15; %convert to kelvin
    %remove -999s through linear interpolation
    temp_var = air_T;
    ind = 1:length(temp_var);
    ix = logical(temp_var > -473.15); %-200 degrees
    temp_var2 = temp_var;
    air_T=interp1(ind(ix),temp_var2(ix),ind,'linear');
    air_T = air_T';

    %SW down
    SW_down = Data{:,'ISWR'};
    SW_down(SW_down < 0) = 0;
    %remove -999s through linear interpolation
    temp_var = SW_down;
    ind = 1:length(temp_var);
    ix = logical(temp_var >= 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    SW_down=interp1(ind(ix),temp_var2(ix),ind,'linear');
    SW_down=SW_down';

    %SW_out
     %SW down
    SW_out = Data{:,'OSWR'};
    SW_out(SW_out < 0) = 0;
    %remove -999s through linear interpolation
    temp_var = SW_out;
    ind = 1:length(temp_var);
    ix = logical(temp_var >= 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    SW_out=interp1(ind(ix),temp_var2(ix),ind,'linear');
    SW_out=SW_out';
    
    %albedo 
%     albedo = SW_out./SW_down;
    % albedo
    albedo = Data{:,'Alb'};
    %remove -999s through linear interpolation
    temp_var = albedo;
    ind = 1:length(temp_var);
    ix = logical(temp_var > 0 & temp_var < 1); %0 to 1 range
    temp_var2 = temp_var;
    albedo=interp1(ind(ix),temp_var2(ix),ind,'linear');
    albedo=albedo';

    %pressure
    pressure = Data{:,'P'};
    pressure(pressure == 0) = NaN;
    pressure = pressure/10; %convert from hPa to kPa
    %remove -999s through linear interpolation
    temp_var = pressure;
    ind = 1:length(temp_var);
    ix = logical(temp_var > 0); %0 chosen as physically impossible
    temp_var2 = temp_var;
    pressure=interp1(ind(ix),temp_var2(ix),ind,'linear');
    pressure=pressure';

    %relative humidity
    rh = Data{:,'RH1'};
    rh(rh == 0) = NaN;
    %remove -999s through linear interpolation
    temp_var = rh;
    ind = 1:length(temp_var);
    ix = logical(temp_var > -473.15); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    rh=interp1(ind(ix),temp_var2(ix),ind,'linear');
    rh=rh';

    %convert relative humidity to specific humidity following https://earthscience.stackexchange.com/questions/5076/how-to-calculate-specific-humidity-with-relative-humidity-temperature-and-pres
    %constants
    Rv = 461.5; %J/kg*K specific gas constant
    Lv = 2257*10^3; %J/kg specific enthalpy of vaporisation
    T0 = 273.15; %K reference temperature, may need changing
    Rd = 287.058; %J/kg*K specific gas constant for dry air
    es0 = 613; %Pa water saturation pressure at 0 degrees

    %rh humidity from % to proportion
    rh_proportion = rh./100;

    %loop
    es = es0.*exp((Lv./Rv)*(1./T0 - 1./air_T));
    e = rh_proportion.*es;
    w = (e.*Rd)./(Rv.*(pressure.*1000 - e)); %*1000 for Kpa to pa
    q = w./(w + 1);

    q(q == 0) = NaN; %get rid of zeros that have crept in
    hum = q; %save specific humidity as hum in kg/kg

    %windspeed
    wind_speed = Data{:,'VW1'};
    wind_speed(wind_speed == 0) = NaN;
    %remove -999s through linear interpolation
    temp_var = wind_speed;
    ind = 1:length(temp_var);
    ix = logical(temp_var > 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    wind_speed=interp1(ind(ix),temp_var2(ix),ind,'linear');
    wind_speed=wind_speed';
    wind_speed(wind_speed<1)=1;

    % 向下长波计算
    %constants
    e_oc = 0.952;           %emissivity of overcast sky
    n = 0.6;                %cloudiness, unitless between 0 and 1
    sigma = 5.67*10^-8;     %stefan boltzmann constant [W/m2-K]
    b = 0.484;              %constant
    m = 8;                  %constant
    p = 4;                  %constant
    e_sat = 610.8.*exp(22.47.*(1 - (273.15./air_T))); %for ice
    %calculate e_air (reference vapour pressure)
    e_air = e_sat.*rh_proportion; %/100 to go from % to proportion
    %calculate e_cs (clear sky emissivity)
    e_cs = 0.23 + b.*((e_air./air_T).^(1./m));
    %calculate LW_in
    LW_down = (e_cs.*(1 - n.^p) + e_oc.*(n.^p)).*sigma.*air_T.^4;
    LW_down(LW_down < 0) = 0;
    %remove -999s through linear interpolation
    temp_var = LW_down;
    ind = 1:length(temp_var);
    ix = logical(temp_var >= 0); %-200 degrees chosen but may need to change, data missing = -999
    temp_var2 = temp_var;
    LW_down=interp1(ind(ix),temp_var2(ix),ind,'linear');
    LW_down=LW_down';

    SW_down_out = SW_down(clip_hour1:clip_hour2);
    SW_down_out(SW_down_out < 0) = 0;
    LW_down_out = LW_down(clip_hour1:clip_hour2);
    LW_down_out(LW_down_out < 0) = 0;
    albedo_out = albedo(clip_hour1:clip_hour2);
    albedo_out(isnan(albedo_out)) = 0;
    air_T_out = air_T(clip_hour1:clip_hour2);
    rh_out = rh(clip_hour1:clip_hour2);
    rh_out(rh_out>100)=100;
    hum_out = hum(clip_hour1:clip_hour2);
    pressure_out = pressure(clip_hour1:clip_hour2);
    wind_speed_out = wind_speed(clip_hour1:clip_hour2);
    time_out = time(clip_hour1:clip_hour2);
end
% spin-up data
clip = yeardays(run_start(1))*24;
SW_down_s = SW_down_out(1:clip);
albedo_s = albedo_out(1:clip);
LW_down_s = LW_down_out(1:clip);
air_T_s = air_T_out(1:clip);
rh_s = rh_out(1:clip);
hum_s = hum_out(1:clip);
pressure_s = pressure_out(1:clip);
wind_speed_s = wind_speed_out(1:clip);
s_data = [SW_down_s,albedo_s,LW_down_s,air_T_s,rh_s,hum_s,pressure_s,wind_speed_s];
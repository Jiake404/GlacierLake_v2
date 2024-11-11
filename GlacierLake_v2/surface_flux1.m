function [q,sens_flux , lat_flux , SW_flux1 , LW_down , LW_out, CT,q1] = surface_flux1(SW_down, albedo_out, air_T, wind_speed, pressure,...
    LW_down, hum, surface_T, albedo, AWS_albedo, Io, e_ice,tur_num)
%GlacierLake subfunction to calculate surface energy flux. Based on Buzzard
%(2017)
        
    if AWS_albedo == 1 %use AWS albedo input
        SW_flux = SW_down*(1 - albedo_out); %by definition of albedo SW_flux = incident light not reflected
    else
        SW_flux = SW_down*(1 - albedo); %use set albedo            
    end

    %only carry forward the SW flux that does not penetrate further into
    %the ice under Beer Lambert law as this is dealt with in the
    %subfunction SW_propagate
    SW_flux1 = SW_flux;
    SW_flux = SW_flux*(1-Io);
    
    %calculate sensible and latent heat following Ebert & Curry (1993) and Buzzard (2017)
    %constants
    g = 9.81;               
    b = 20;                 
    dz = 2;        
    CT0 = 1.3*10^(-3);      
    c = 50.986;     %wubbold, this is different to the value of c in Sammie's code
    air_ro = 1.275; %(kgm^-3) density of dry air
    air_cp = 1005;  %specific heat capacity of dry air
    

    %Richardson number
    Ri = (g*(air_T - surface_T)*dz)/(air_T*wind_speed.^2);
     %CT (function of atmospheric stability)
        if Ri < 0
            CT = CT0*(1 - (2*b*Ri)/(1 + c*abs(Ri.^0.5)));
        else
            CT = CT0*(1 + b*Ri).^(-2);
        end
    if tur_num==3
        Lv = 2.5*10^6;%(Jkg^-1) latent heat of vaporisation
    else
        Lv = 2.83*10^6;%(Jkg^-1) latent heat of vaporisation
    end
    %sensible and latent heat equations
    p_v = 2.53*10.^8.*exp(-5420./surface_T); %wubbold not discussed in Buzzard (2017)
    surface_hum = (0.622.*p_v)./(pressure - 0.378.*p_v); %wubbold not discussed in Buzzard (2017)
    sens_flux = air_ro.*air_cp.*CT*wind_speed.*(air_T - surface_T);
    lat_flux = air_ro.*Lv.*CT.*wind_speed.*(hum - surface_hum); %hum is air humidity  

    %longwave up 
    sigma = 5.67*10^-8; %(W/m2/K) Stefan Boltzmann constant
    LW_out = e_ice*sigma*(surface_T^4);

    q = sens_flux + lat_flux + SW_flux + LW_down - LW_out; %(J/s) surface energy absorption
    q1 = sens_flux + lat_flux + SW_flux1 + LW_down - LW_out; %(J/s) surface energy absorption


end
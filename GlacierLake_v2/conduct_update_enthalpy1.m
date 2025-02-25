function [enthalpy_out, lake_in_e_out] = conduct_update_enthalpy1(k, c, T,...
    grid_profile, enthalpy, t_step, total_grid_num, ro_water, lower_boundary,...
    model_stage, snow_T, snow_k, snow_c, snow_z, snow_e)
%GlacierLake subfunction for heat conduction equation and conversion back
%from temperature to enthalpy 
    
    if model_stage == 1
        
        %calculate thermal diffusivity from thermal conductivity
        alpha = (k./c)/ro_water;

        %initialise vectors for tridiagonal matrix, see paper for equations
        alpha_av = (alpha(1:end - 1) + alpha(2:end))./2; %average thermal conductivity values
        S = (alpha_av.*t_step*3600)./((0.5.*grid_profile(1:end - 1) + 0.5.*grid_profile(2:end)).^2); %*3600 is to go from hours to seconds
        diag_bot = zeros(total_grid_num - 1,1); %bottom diagonal
        diag_bot(1:end - 1) = -S(1:end - 1);
        diag_bot(end) = -2*S(end); %bottom neumann
        diag_mid = zeros(total_grid_num,1); %middle diagonal
        diag_mid(1) = 1; %top dirichlet boundary
        diag_mid(2:end - 1) = 1 + S(1:end - 1) + S(2:end);
        diag_mid(end) = 1 + 2*S(end); %bottom neumann boundary
        diag_top = zeros(total_grid_num - 1,1); %top diagonal
        diag_top(1) = 0; %top dirichlet boundary
        diag_top(2:end) = -S(2:end);

        %create tridiagonal matrix
        A = diag(diag_bot, -1) + diag(diag_mid) + diag(diag_top, 1);

        %known array
        T_n1 = zeros(total_grid_num,1);
        T_n1(1:end - 1) = T(1:end - 1);
        T_n1(end) = lower_boundary; %dirichlet boundary

        %solve matrix equation
        T_out = A\T_n1;

        %convert back to enthalpy for model use
        if model_stage == 1 %entire profile is ice
            enthalpy_out = T_out.*grid_profile.*c.*ro_water; %enthalpy where temperature is below freezing
        else
            d_q = (T_out - T).*c.*grid_profile.*ro_water;
            enthalpy_out = enthalpy + d_q;   
        end
    
    elseif model_stage == 2
        
        %if snow is present concantenate, otherwise don't
        if snow_z > 0
            %concantenate temporary arrays
            T_temp = [snow_T; T];
            k_temp = [snow_k; k];
            c_temp = [snow_c; c];
            z_temp = [snow_z; grid_profile];
            enthalpy_temp = [snow_e; enthalpy];
            
            b = 1; %add in factor for tridiag arrays
        else
            T_temp = T;
            k_temp = k;
            c_temp = c;
            z_temp = grid_profile;
            enthalpy_temp = enthalpy;    
            
            b = 0; %add in factor for tridiag arrays
        end
            
        
        %calculate thermal diffusivity from thermal conducvity 
        alpha = (k_temp./c_temp)/ro_water;
            
        %initialise vectors for tridiagonal matrix, see thesis for equations
        alpha_av = (alpha(1:end - 1) + alpha(2:end))./2; %average thermal conductivity values
        S = (alpha_av.*t_step*3600)./((0.5.*z_temp(1:end - 1) + 0.5.*z_temp(2:end)).^2); %*3600 is to go from hours to seconds
        diag_bot = zeros(total_grid_num - 1 + b,1); %bottom diagonal
        diag_bot(1:end - 1) = -S(1:end - 1);
        diag_bot(end) = -2*S(end); %bottom neumann
        diag_mid = zeros(total_grid_num + b,1); %middle diagonal
        diag_mid(1) = 1; %top dirichlet boundary
        diag_mid(2:end - 1) = 1 + S(1:end - 1) + S(2:end);
        diag_mid(end) = 1 + 2*S(end); %bottom neumann boundary
        diag_top = zeros(total_grid_num - 1 + b,1); %top diagonal
        diag_top(1) = 0; %top dirichlet boundary
        diag_top(2:end) = -S(2:end);

        %create tridiagonal matrix
        A = diag(diag_bot, -1) + diag(diag_mid) + diag(diag_top, 1);

        %known array
        T_n1 = zeros(total_grid_num + b,1);
        T_n1(1:end - 1) = T_temp(1:end - 1);
        T_n1(end) = lower_boundary; %dirichlet boundary

        %solve matrix equation
        T_out = A\T_n1;
        
        %change in enthalpy 
        d_q = (T_out - T_temp).*c_temp.*z_temp.*ro_water;
        enthalpy_out = enthalpy_temp + d_q;      
        
        lake_in_e_out = 0; %matlab is a hungry beast
        
    end
    
    
end


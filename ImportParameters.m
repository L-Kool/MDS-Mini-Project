function params = ImportParameters ()
    %% Defining all parameters

    % General parameters
    params.rho_w = 1000;
    params.cp = 4186;
    params.g = 9.81;
    params.T_env = 280.15; % 17 C
    
    % LUMI parameters
    params.V_w_lumi = 10;          
    params.C_lumi = params.rho_w * params.cp * params.V_w_lumi;
    params.R_ex_lumi = 2.8e-4;
    
    % Building parameters 
    params.V_water_B = 3;
    params.C_B_office = 5e7;
    params.C_B_res = 2e7;
    params.R_B_office_0 = 0.006;
    params.R_B_res_0 = 0.008;
    params.C_rad = params.rho_w * params.cp * params.V_water_B; 
    params.R_rad_low = 0.008;  
    params.R_rad_high = 0.5;
    params.R_ex_B_low = 0.0005;
    params.R_ex_B_high = 0.2;
    params.R_wall_ij = 0.02;
    
    % Pump & Circuit parameters
    params.R_dh_flow = 1e7; %1e7;
    params.k_pump = 2e-4; %6e-4 1e-3
    params.K_dh_motor = 3.8;
    params.R_dh_motor = 0.6;
    params.C_dh_motor = 1e-3;
    params.eta_pump = 0.75;
    params.k_load = (params.R_dh_flow * params.k_pump^2) / ...
                    (params.eta_pump * params.K_dh_motor^2);
    
    % Hydro & Battery parameters
    params.A_res = 3000;
    params.R_hydro_flow = 1.2e5; %7
    params.h_hydro_min = 10;
    params.h_hydro_max = 25;
    params.eta_hydro_turbine = 0.85;
    params.C_hydro_battery = 1.08e6;
    params.V_batt = 600;
    params.P_grid = 15000;
    params.I_charging = params.P_grid / params.V_batt;
    
    % Load Piping parameters
    params.L_dh = 1500;
    params.A = 0.08;
    params.N = 100;
    params.Delta_x = params.L_dh / params.N;
    params.k_dh = 0.4;
    params.G_dh_thermal = 0.3; 
    params.C_seg = params.rho_w * params.A * params.cp;
    params.k_prime = params.k_dh / (params.rho_w * params.cp);
    params.alpha = params.G_dh_thermal / (params.C_seg);

end 
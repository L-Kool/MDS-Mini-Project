% Main script to simulate the complete system

% Importing all parameters
params = ImportParameters();

%% Defining simulation conditions

% Defining time span
t_span = [0 3600*48]; % 48 hours

% Imposing that pumps are constantly ON
params.SIM.Pump_B_State = 'ON';
params.SIM.Pump_Sp_State = 'CLOSED';

% Providing constant flow to reservoir so it will not deplete
w_in_sim = [1.8 1.2];
params.SIM.w_in = w_in_sim; % Constant river flow

Q_lumi_sim = [1.9e5 1.2e5];
params.SIM.Q_lumi = Q_lumi_sim; % 2MW constant heat load

% Defining locations of building along segmented pipe
params.pipe_nodes.B1 = params.N / 5;
params.pipe_nodes.B2 = 2 * params.N/5; 
params.pipe_nodes.B3 = 3 * params.N/5;
params.pipe_nodes.B4 = 4 * params.N/5;

% Defining initial state of super vector
x0 = zeros(12 + params.N, 1);
x0(1) = 353.15; % T_LUMI (80 C) 353
x0(2:2:8) = 293.15; % T_B (20 C)
x0(3:2:9) = 303.15; % T_rad (30 C)
x0(10) = 590;     % Vm (steady state)
x0(11) = 12;      % h_res (12m)
x0(12) = params.C_hydro_battery; % q_batt (full)
x0(13:end) = 318.15; % T_pipe (70 C)

%% Running the simulations
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'MaxStep', 10, 'Stats', 'on');
[t_sol, x_sol] = ode15s(@(t, x) complete_dynamics(t, x, params), ...
                                    t_span, x0, options);

%% Post-process, calculating outputs
h_sol = x_sol(:, 11);
w_out_sol = zeros(size(t_sol));
P_mech = zeros(size(t_sol)); % Generator Power

for i = 1:length(t_sol)
    h_i = h_sol(i);
    if h_i > params.h_hydro_min
        P_fluid = (params.rho_w * params.g * h_i)^2 / params.R_hydro_flow;
        w_out_sol(i) = P_fluid / (params.rho_w * params.g * h_i);
        % Calculate mechanical power
        P_mech(i) = params.eta_hydro_turbine * P_fluid;
    end
end

Vm_sol = x_sol(:, 10);       

% Calculate omega_pump at each time step
omega_pump_sol = Vm_sol / params.K_dh_motor;

% Calculate flow rate w_dh at each time step [cite: 150]
w_dh_sol = params.k_pump * omega_pump_sol;

% Calculate velocity v at each time step [cite: 75]
v_sol = w_dh_sol / params.A;


%% Plotting results
close all;
figure(1);
subplot(2,1,1);
plot(t_sol/3600, x_sol(:, 1) - 273.15, 'LineWidth', 3);
title('LUMI Temperature', 'FontSize', 14); 
xlabel('Time (hr)', 'FontSize', 14); 
ylabel('Temp (C)', 'FontSize', 14);
xlim([0 48]);
grid on;

% Pipe inlet and outlet temperature
subplot(2,1,2);
T_pipe_in = x_sol(:,13) - 273.15;
T_pipe_out = x_sol(:, end) - 273.15;
plot(t_sol/3600, T_pipe_out, 'LineWidth', 3);
hold on;
plot(t_sol/3600, T_pipe_in, 'LineWidth', 3);
title('Pipe Inlet and Outlet Temp', 'FontSize', 14);
xlabel('Time (hr)', 'FontSize', 14); 
ylabel('Temp (C)', 'FontSize', 14);
xlim([0 48]);
legend('Outlet', 'Inlet', 'FontSize', 14)
grid on;

% Building temperature
figure(2);
subplot(2,1,1);
plot(t_sol/3600, x_sol(:, 2) - 273.15, 'r', 'DisplayName', 'B1 (Office)', 'LineWidth', 3);
hold on;
plot(t_sol/3600, x_sol(:, 4) - 273.15, 'b', 'DisplayName', 'B2 (Res)', 'LineWidth', 3);
plot(t_sol/3600, x_sol(:, 6) - 273.15, 'g', 'DisplayName', 'B3 (Office)', 'LineWidth', 3);
plot(t_sol/3600, x_sol(:, 8) - 273.15, 'm', 'DisplayName', 'B4 (Res)', 'LineWidth', 3);
title('Building Temperatures', 'FontSize', 14);
xlabel('Time (hr)', 'FontSize', 14); 
ylabel('Temp (C)', 'FontSize', 14);
xlim([0 48]);
legend('FontSize', 14);
grid on;

% Radiator temperature
subplot(2,1,2);
T_rad_1 = x_sol(:,3) - 273.15;
T_rad_2 = x_sol(:,5) - 273.15;
T_rad_3 = x_sol(:,7) - 273.15;
T_rad_4 = x_sol(:,9) - 273.15;
plot(t_sol/3600, T_rad_1, 'r', 'DisplayName', 'Radiator B1 (Office', 'LineWidth', 3);
hold on
plot(t_sol/3600, T_rad_2, 'b', 'DisplayName', 'Radiator B2 (Res)', 'LineWidth', 3);
plot(t_sol/3600, T_rad_3, 'g', 'DisplayName', 'Radiator B3 (Office)', 'LineWidth', 3);
plot(t_sol/3600, T_rad_4, 'm', 'DisplayName', 'Radiator B4 (Res)', 'LineWidth', 3);
title('Radiator Temperatures', 'FontSize', 14);
xlabel('Time (hr)', 'FontSize', 14); 
ylabel('Temp (C)', 'FontSize', 14);
xlim([0 48]);
legend('FontSize', 14)
grid on;

% Plotting reservoir height and mech. power
figure(3)
subplot(2,1,1);
plot(t_sol/3600, x_sol(:, 11), 'DisplayName', 'Reservoir height', 'LineWidth', 3);
title('Reservoir Height', 'FontSize', 14);
xlabel('Time (hr)', 'FontSize', 14); 
ylabel('Height (m)', 'FontSize', 14);
xlim([0 48]);
legend('FontSize', 14)
grid on;

subplot(2,1,2);
plot(t_sol/3600, P_mech, 'DisplayName', 'Turbine mechanical power', 'LineWidth', 3);
title('Turbine mechanical power', 'FontSize', 14);
xlabel('Time (hr)', 'FontSize', 14); 
ylabel('Power (W)', 'FontSize', 14);
xlim([0 48]);
legend('FontSize', 14)
legend; grid on;

% Plotting temperature profiles
% Extracing pipe temperatures
T_pipe_sol = x_sol(:, 13:end); 
times_sec = t_sol;            
num_segments = params.N;
pipe_length = params.L_dh;
delta_x = params.Delta_x;

% Determinig centres of pipe segments
x_pipe = linspace(delta_x/2, pipe_length - delta_x/2, num_segments); 

% Selecting time points to plot
num_time_points = 5;
time_indices = round(linspace(1, length(times_sec), num_time_points));
selected_times_hr = times_sec(time_indices) / 3600; 

% Plotting
figure(4);
hold on; 

for i = 2:num_time_points
    idx = time_indices(i);
    temp_profile_C = T_pipe_sol(idx, :) - 273.15; 
    plot(x_pipe, temp_profile_C, 'LineWidth', 2, ...
         'DisplayName', sprintf('t = %.1f hr', selected_times_hr(i)));
end

hold off;
title('Pipe Temperature Profile Evolution', 'FontSize', 14);
xlabel('Pipe Position (m)', 'FontSize', 14);
ylabel('Temperature (°C)', 'FontSize', 14);
legend('FontSize', 14);
grid on;

% Adding building locations
y_limits = ylim; 
building_nodes = [params.pipe_nodes.B1, params.pipe_nodes.B2, params.pipe_nodes.B3, params.pipe_nodes.B4];
building_pos = x_pipe(building_nodes); 

for i = 1:length(building_pos)
    line([building_pos(i) building_pos(i)], y_limits, 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off', 'LineWidth', 2);
    text(building_pos(i), y_limits(1) + 0.05*diff(y_limits), sprintf(' B%d', i), 'HorizontalAlignment', 'left');
end
ylim(y_limits); 

%% =================================================================
%  Function for complete model
% x_total = [
%   x(1)       % T_LUMI      
%
%   x(2)       % T_B1        
%   x(3)       % T_rad1      
%
%   x(4)       % T_B2        
%   x(5)       % T_rad2      
%
%   x(6)       % T_B3       
%   x(7)       % T_rad3      
%
%   x(8)       % T_B4       
%   x(9)       % T_rad4      
%
%   x(10)      % Vm          
%   x(11)      % h_res       
%   x(12)      % q_batt      
%
%   x(13)      % T_pipe(1)   
%   ...        % ...
%   x(12+N)    % T_pipe(N)  
% ]
function x_dot = complete_dynamics(t, x_total, params)
    
    % Unpacking state vector
    T_lumi = x_total(1);
    
    T_B1 = x_total(2);
    T_rad1 = x_total(3);
    
    T_B2 = x_total(4);
    T_rad2 = x_total(5);
    
    T_B3 = x_total(6);
    T_rad3 = x_total(7);
    
    T_B4 = x_total(8);
    T_rad4 = x_total(9);
    
    Vm = x_total(10);
    h_res = x_total(11);
    q_batt = x_total(12);
    
    T_pipe = x_total(13 : 12+params.N);
    
    % Calculating coupling equations
    
    % Pump & pipe coupling 
    omega_pump = Vm / params.K_dh_motor;
    w_dh = params.k_pump * omega_pump;
    
    % Lumi & pipe coupling 
    T_dh_return = T_pipe(end);
    Q_ex_lumi = (T_lumi - T_dh_return) / params.R_ex_lumi;
    T_dh_inlet = T_dh_return + Q_ex_lumi / (w_dh * params.rho_w * params.cp + 1e-3);

    % Building & pipe coupling (heat sinks) 
    T_pipe_B1 = T_pipe(params.pipe_nodes.B1);
    T_pipe_B2 = T_pipe(params.pipe_nodes.B2);
    T_pipe_B3 = T_pipe(params.pipe_nodes.B3);
    T_pipe_B4 = T_pipe(params.pipe_nodes.B4);
    
    % (Using R_ex_B_low because pump state is ON)
    Q_sink_B1 = (T_pipe_B1 - T_rad1) / params.R_ex_B_low; 
    Q_sink_B2 = (T_pipe_B2 - T_rad2) / params.R_ex_B_low;
    Q_sink_B3 = (T_pipe_B3 - T_rad3) / params.R_ex_B_low;
    Q_sink_B4 = (T_pipe_B4 - T_rad4) / params.R_ex_B_low;
    
    Q_sinks_vec = zeros(params.N, 1);
    Q_sinks_vec(params.pipe_nodes.B1) = Q_sink_B1;
    Q_sinks_vec(params.pipe_nodes.B2) = Q_sink_B2;
    Q_sinks_vec(params.pipe_nodes.B3) = Q_sink_B3;
    Q_sinks_vec(params.pipe_nodes.B4) = Q_sink_B4;
    
    % Supervisor logic 
    pump_B_state = params.SIM.Pump_B_State;
    pump_Sp_state = params.SIM.Pump_Sp_State;
    
    if h_res > params.h_hydro_min
        V_in = 600; 
        batt_mode = 'IDLE';
    else
        V_in = 600; 
        batt_mode = 'DISCHARGE'; 
    end
    
    % Computing dynamics
    
    % LUMI
    inputs_lumi.Q_lumi = params.SIM.Q_lumi;
    inputs_lumi.T_return = T_dh_return; 
    T_lumi_dot = SystemDynamics.lumi_dynamics(t, T_lumi, params, inputs_lumi);
    
    % Building 1
    inputs_B1.T_env = params.T_env;
    inputs_B1.T_neighbor = T_B2; 
    inputs_B1.T_dh = T_pipe_B1;
    x_dot_B1 = SystemDynamics.building_dynamics(t, [T_B1; T_rad1], params, inputs_B1, 'Office', pump_B_state); 
    
    % Building 2
    inputs_B2.T_env = params.T_env;
    inputs_B2.T_neighbor = T_B1; 
    inputs_B2.T_dh = T_pipe_B2;
    x_dot_B2 = SystemDynamics.building_dynamics(t, [T_B2; T_rad2], params, inputs_B2, 'Residential', pump_B_state); 
    
    % Building 3
    inputs_B3.T_env = params.T_env;
    inputs_B3.T_neighbor = T_B4; 
    inputs_B3.T_dh = T_pipe_B3;
    x_dot_B3 = SystemDynamics.building_dynamics(t, [T_B3; T_rad3], params, inputs_B3, 'Office', pump_B_state); 
    
    % Building 4
    inputs_B4.T_env = params.T_env;
    inputs_B4.T_neighbor = T_B3;
    inputs_B4.T_dh = T_pipe_B4;
    x_dot_B4 = SystemDynamics.building_dynamics(t, [T_B4; T_rad4], params, inputs_B4, 'Residential', pump_B_state); 
    
    % Pump Circuit
    inputs_pump.V_in = V_in;
    inputs_pump.Sp_state = pump_Sp_state;
    Vm_dot = SystemDynamics.pump_circuit_dynamics(t, Vm, params, inputs_pump);
    
    % Reservoir
    inputs_res.w_in = params.SIM.w_in;
    h_res_dot = SystemDynamics.reservoir_dynamics(t, h_res, params, inputs_res);
    
    % Battery
    inputs_batt.Vm = Vm;
    inputs_batt.mode = batt_mode;
    q_batt_dot = SystemDynamics.battery_dynamics(t, q_batt, params, inputs_batt);

    % Pipe
    inputs_pipe.v = w_dh / params.A;
    inputs_pipe.T_inlet = T_dh_inlet;
    inputs_pipe.T_env = params.T_env;
    inputs_pipe.Q_sinks = Q_sinks_vec;
    T_pipe_dot = SystemDynamics.piping_dynamics(t, T_pipe, params, inputs_pipe);
    
    % Packing the state vector
    x_dot = [
        T_lumi_dot;
        x_dot_B1;
        x_dot_B2;
        x_dot_B3;
        x_dot_B4;
        Vm_dot;
        h_res_dot;
        q_batt_dot;
        T_pipe_dot
    ];
end
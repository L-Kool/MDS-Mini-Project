% Main script to simulate a single building
clear; clc; close all;

% Importing parameters
params = ImportParameters();

% Simulation conditions
t_span = [0 3600*48]; % Simulate for 48 hours
x0 = [293.15; 303.15]; % Initial T_B = 20C, T_rad = 30C

% Defining inputs
inputs.T_env = 288.15;      % 15 C
inputs.T_neighbor = 288.15; % 25 C (neighboring building)
inputs.T_dh = 343.15;       % 70 C (district heating supply)

% Run simulation for office building, pumps ON
[t_sol_office_on, x_sol_office_on] = ...
    ode45(@(t, x) SystemDynamics.building_dynamics(t, x, params, inputs, 'Office', 'ON'), t_span, x0);

% Run simulation for residential building, pumps ON
[t_sol_res_on, x_sol_res_on] = ...
    ode45(@(t, x) SystemDynamics.building_dynamics(t, x, params, inputs, 'Residential', 'ON'), t_span, x0);

% Run simulation for office building, pumps OFF
[t_sol_office_off, x_sol_office_off] = ...
    ode45(@(t, x) SystemDynamics.building_dynamics(t, x, params, inputs, 'Office', 'OFF'), t_span, x0);

% Run simulation for residential building, pumps OFF
[t_sol_res_off, x_sol_res_off] = ...
    ode45(@(t, x) SystemDynamics.building_dynamics(t, x, params, inputs, 'Residential', 'OFF'), t_span, x0);

% Plotting results (Office, pumps ON)
figure;
subplot(2,2,1);
plot(t_sol_office_on / 3600, x_sol_office_on - 273.15, 'LineWidth', 2);
title('Building 1 (Office), pump ON ', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Temperature (C)', 'FontSize', 14);
xlim([0 48]);
legend('T_{Building}', 'T_{Radiator}', 'FontSize', 14);
grid on;

% Plotting results (Residential, pumps ON)
subplot(2,2,2);
plot(t_sol_res_on / 3600, x_sol_res_on - 273.15, 'LineWidth', 2);
title('Building 2 (Residential), pump ON ', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Temperature (C)', 'FontSize', 14);
xlim([0 48]);
legend('T_{Building}', 'T_{Radiator}', 'FontSize', 14);
grid on;

% Plotting results (Office, pumps OFF)
subplot(2,2,3);
plot(t_sol_office_off / 3600, x_sol_office_off - 273.15, 'LineWidth', 2);
title('Building 1 (Office), pumps OFF ', 'FontSize', 14);
xlabel('Time (hours)','FontSize', 14);
ylabel('Temperature (C)', 'FontSize', 14);
xlim([0 48]);
legend('T_{Building}', 'T_{Radiator}', 'FontSize', 14);
grid on;

% Plotting results (Residential, pumps OFF)
subplot(2,2,4);
plot(t_sol_res_off / 3600, x_sol_res_off - 273.15, 'LineWidth', 2);
title('Building 2 (Residential), pumps OFF ', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Temperature (C)', 'FontSize', 14);
xlim([0 48]);
legend('T_{Building}', 'T_{Radiator}', 'FontSize', 14);
grid on;

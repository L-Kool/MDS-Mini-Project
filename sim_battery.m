% Main script to simulate the backup battery
clear; clc; close all;

% Importing parameters
params = ImportParameters();

% Simulation 1: Discharging
t_span_dis = [0 3600 * 12]; % 12 hours
x0_dis = [params.C_hydro_battery]; 

% Inputs for discharge test
inputs_dis.Vm = 580; 
inputs_dis.mode = 'DISCHARGE';

[t_sol_dis, x_sol_dis] = ode45(@(t, x) SystemDynamics.battery_dynamics(t, x, params, inputs_dis), ...
                                    t_span_dis, x0_dis);

% Simulation 2: Charging
t_span_ch = [0 3600 * 12]; % 12 hours
x0_ch = [x_sol_dis(end)];  % Start from the discharged state

% Inputs for charge test
inputs_ch.Vm = 0; % Vm is irrelevant during charging
inputs_ch.mode = 'CHARGE';

[t_sol_ch, x_sol_ch] = ode45(@(t, x) SystemDynamics.battery_dynamics(t, x, params, inputs_ch), ...
                                    t_span_ch, x0_ch);

% Plotting results
% Convert charge (Coulombs) to State of Charge (%)
SoC_dis = 100 * x_sol_dis / params.C_hydro_battery;
SoC_ch = 100 * x_sol_ch / params.C_hydro_battery;

% Combine results for one plot
t_final = [t_sol_dis; t_sol_ch + t_sol_dis(end)];
SoC_final = [SoC_dis; SoC_ch];

close all;
figure;
plot(t_final / 3600, SoC_final, 'LineWidth', 2);
title('Battery State of Charge (SoC)', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('SoC (%)', 'FontSize', 14);
grid on;
ylim([-1 110]);
xlim([0 24]);
xline(t_sol_dis(end)/3600, 'r--', 'Label', 'Start Charging', 'LineWidth', 3);
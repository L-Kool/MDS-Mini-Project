% Main script to simulate the LUMI component in isolation.
clear; clc; close all;

% Importing parameters
params = ImportParameters();

% Defining simulation conditions
t_span = [0 3600 *48]; % Simulate for 1 hour
T_initial = 303.15; % Initial temp (30 C)
x0 = T_initial;

% Defining inputs for isolated simulation
inputs.Q_lumi = [1e5 0.8e5];  
inputs.T_return = 313.15; % 40 C 

% Running the simulation
[t_sol, x_sol] = ode45(@(t, x) SystemDynamics.lumi_dynamics(t, x, params, inputs), t_span, x0);
% Plotting the results
close all;
figure;
plot(t_sol / 3660, x_sol(:, 1) - 273.15, 'LineWidth', 3);
title('LUMI Cooling Circuit Temperature', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
xlim([0 48]);
ylim([20 75]);
ylabel('Temperature T_{LUMI} (C)', 'FontSize', 14);
legend('LUMI Temperature', 'FontSize', 14);
grid on;
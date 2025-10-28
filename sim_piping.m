% Main script to simulate the PDE piping model
clear; clc; close all;

% Importing parameters
params = ImportParameters();

% Simulation and initial conditions
t_span = [0 3600*24]; 
x0 = ones(params.N, 1) * params.T_env; 

% Defining inputs
w_dh_input = 0.03; % Flow of water in district piping
inputs.v = w_dh_input / params.A; 
inputs.T_inlet = 343.15; % 70 C
inputs.T_env = params.T_env;

% Defining sinks for isolated tests
nodeB1 = params.N / 5;
nodeB2 = 2 * nodeB1;
nodeB3 = 3 * nodeB1;
nodeB4 = 4 * nodeB1;

% Initialize the heat sink values for the nodes
inputs.Q_sinks = zeros(params.N, 1); % No buildings attached
inputs.Q_sinks(nodeB1) = 14000; % Heat sink at nodeB1
inputs.Q_sinks(nodeB2) = 5000; % Heat sink at nodeB2
inputs.Q_sinks(nodeB3) = inputs.Q_sinks(nodeB1); % Heat sink at nodeB3
inputs.Q_sinks(nodeB4) = inputs.Q_sinks(nodeB2); % Heat sink at nodeB4

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'MaxStep', 7, 'Stats', 'on');

% Running simulation
[t_sol, x_sol] = ode15s(@(t, x) SystemDynamics.piping_dynamics(t, x, params, inputs), ...
                                    t_span, x0, options);

% Plotting results
figure;
x_pipe = linspace(params.Delta_x, params.L_dh, params.N);
plot(x_pipe, x_sol(end, :) - 273.15, 'LineWidth', 3);
title('Piping Temperature Profile at t = 24 hours', 'FontSize', 14);
xlabel('Pipe Position (m)', 'FontSize', 14);
ylabel('Temperature (C)', 'FontSize', 14);
grid on;
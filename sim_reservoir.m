% Main script to simulate the hydropower reservoir
clear; clc; close all;

% Importing parameters
params = ImportParameters();

% Defining simulation conditions
t_span = [0 3600*48]; 
x0 = 12;               

% Defining inputs
t_sim = linspace(t_span(1), t_span(end), 1000);
inputs.w_in = [1.8 1.2];


% Run simulation
[t_sol, x_sol] = ...
    ode15s(@(t, x) SystemDynamics.reservoir_dynamics(t, x, params, inputs), t_span, x0);

% Post-process, calculate outputs
h_sol = x_sol(:, 1);
w_out_sol = zeros(size(t_sol));
P_mech = zeros(size(t_sol)); % Generator Power

for i = 1:length(t_sol)
    h_i = h_sol(i);
    if h_i > params.h_hydro_min
        P_fluid = (params.rho_w * params.g * h_i)^2 / params.R_hydro_flow;
        w_out_sol(i) = P_fluid / (params.rho_w * params.g * h_i);
        
        % Calculate mechanical power
        P_mech(i) = params.eta_hydro_turbine * P_fluid;
    else
        w_out_sol(i) = 0; 
        P_mech(i) = 0;    
    end
end

% Plotting results
close all;
figure;
subplot(3, 1, 1); 
plot(t_sol / 3600, h_sol, 'b-', 'LineWidth', 2);
hold on;
yline(params.h_hydro_min, 'r--', 'LineWidth', 2);
yline(params.h_hydro_max, 'r--', 'LineWidth', 2);
title('Reservoir Level h(t)', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Height (m)', 'FontSize', 14);
ylim([9 26])
xlim([0 48]);
legend('Water Level', 'Min/Max Levels', 'FontSize', 14);
grid on;

subplot(3, 1, 2);
w_in_1 = inputs.w_in(1)*ones(ceil(size(t_sol)/2));
w_in_2 = inputs.w_in(2)*ones(round(size(t_sol)/2));
w_in_array = horzcat(w_in_1, w_in_2);
plot(t_sol / 3600, w_in_array(1:end-1), 'DisplayName', 'Inflow (w_{in})', 'LineWidth', 2);
hold on;
plot(t_sol / 3600, w_out_sol, 'DisplayName', 'Outflow (w_{out})', 'LineWidth', 2);
ylim([0.5 2]);
title('Reservoir Flow Rates', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Flow Rate (m^3/s)', 'FontSize', 14);
legend('FontSize', 14);
grid on;

subplot(3, 1, 3); 
plot(t_sol / 3600, P_mech / 1000, 'LineWidth', 2);
title('Generator Mechanical Power Output', 'FontSize', 14);
xlabel('Time (hours)', 'FontSize', 14);
ylabel('Power (kW)', 'FontSize', 14);
grid on;
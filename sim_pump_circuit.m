% Main script to simulate the pump electrical circuit
clear; clc; close all;

% Importing parameters
params = ImportParameters();

% Parameters for simulation (pump ON)
t_span_on = [0 3600]; 
x0_on = 0;       
inputs_on.V_in = 600;         
inputs_on.Sp_state = 'CLOSED'; 

% Running simulation with pumps ON
[t_sol_on, x_sol_on] = ode15s(@(t, x) SystemDynamics.pump_circuit_dynamics(t, x, params, inputs_on), ...
                                    t_span_on, x0_on);

% Parameters for simulation (pump OFF, spin-down)
t_span_off = [0 3600];
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'Stats', 'on');
x0_off = [x_sol_on(end)]; 

inputs_off.V_in = 0;         
inputs_off.Sp_state = 'OPEN';

% Running simulation with pump OFF
[t_sol_off, x_sol_off] = ode45(@(t, x) SystemDynamics.pump_circuit_dynamics(t, x, params, inputs_off), ...
                                    t_span_off, x0_off);

% Plotting results
close all;
figure;
subplot(2, 1, 1);
plot(t_sol_on/3600, x_sol_on(:, 1), 'b-', 'LineWidth', 2);
hold on;
plot((t_sol_off + t_span_on(end))/3600, x_sol_off(:, 1), 'r--', 'LineWidth', 2);
title('Motor Voltage V_m', 'FontSize', 14);
xlabel('Time (h)', 'FontSize', 14);
ylabel('Voltage (V)', 'FontSize', 14);
ylim([-1 620]);
legend('Turn-On (Sp=1)', 'Turn-Off (Sp=0)', 'FontSize', 14);
grid on;

omega_on = x_sol_on(:, 1) / params.K_dh_motor;
omega_off = x_sol_off(:, 1) / params.K_dh_motor;

subplot(2, 1, 2);
plot(t_sol_on/3600, omega_on, 'b-', 'LineWidth', 2);
hold on;
plot((t_sol_off + t_span_on(end))/3600, omega_off, 'r--', 'LineWidth', 2);
title('Pump Speed \omega_{pump}', 'FontSize', 14);
xlabel('Time (h)', 'FontSize', 14);
ylabel('Speed (rad/s)', 'FontSize', 14);
ylim([-10 170]);
legend('Turn-On (Sp=1)', 'Turn-Off (Sp=0)', 'FontSize', 14);
grid on;
classdef SystemDynamics
% A class created to have all the dynamic models in one place
% They are called in the separate simulation files

    methods (Static)
        %% 1. LUMI COOLING CIRCUIT 
        function T_lumi_dot = lumi_dynamics(t, T_lumi, params, inputs)
            
            % Unpacking parameters
            C_lumi = params.C_lumi;         
            R_ex_lumi = params.R_ex_lumi;
            
            % Unpacking inputs            
            T_return = inputs.T_return; 
            if t <= 24*3600
                Q_lumi = inputs.Q_lumi(1);
            else 
                Q_lumi = inputs.Q_lumi(2);
            end

            % Model
            Q_out = (T_lumi - T_return) / R_ex_lumi;
            T_lumi_dot = (Q_lumi - Q_out) / C_lumi;

        end
        
        %% 2. BUILDING MODEL
        function x_dot = building_dynamics(t, x, params, inputs, building_type, pump_state)
            
            % Unpacking state
            T_Bi = x(1);    
            T_rad_i = x(2); 
            
            % Unpacking parameters
            C_B_office = params.C_B_office;
            C_B_res = params.C_B_res;
            R_B_office_0 = params.R_B_office_0;
            R_B_res_0 = params.R_B_res_0;
            C_rad = params.C_rad;
            R_wall_ij = params.R_wall_ij;
            
            % Unpacking inputs
            T_env = inputs.T_env;
            T_neighbor = inputs.T_neighbor;
            T_dh = inputs.T_dh;          
            
            if strcmpi(pump_state, 'ON') 
                R_rad = params.R_rad_low;
                R_ex_B = params.R_ex_B_low;
            elseif strcmpi(pump_state, 'OFF') 
                R_rad = params.R_rad_high;
                R_ex_B = params.R_ex_B_high;
            end

            if strcmpi(building_type, 'Office')
                C_Bi = C_B_office;
                R_B_i_0 = R_B_office_0;
            elseif strcmpi(building_type, 'Residential')
                C_Bi = C_B_res;
                R_B_i_0 = R_B_res_0;
            end 
                
            % Model of building air
            Q_rad = (T_rad_i - T_Bi) / R_rad;
            Q_out_env = (T_Bi - T_env) / R_B_i_0;
            Q_building = (T_Bi - T_neighbor) / R_wall_ij; 
            
            T_Bi_dot = (Q_rad - Q_out_env - Q_building) / C_Bi;
            
            % Model of radiator water circuit
            Q_dh_in = (T_dh - T_rad_i) / R_ex_B;
            
            T_rad_i_dot = (Q_dh_in - Q_rad) / C_rad;
            
            % Combining models
            x_dot = [T_Bi_dot; T_rad_i_dot];
        end
        
        %% 3. DISTRICT HEATING PIPING 
        function T_dot = piping_dynamics(t, T, params, inputs)

            % Unpack parameters
            N = params.N;
            Delta_x = params.Delta_x;
            k_prime = params.k_prime;
            alpha = params.alpha;
            C_seg = params.C_seg;
        
            % Unpack inputs
            v = inputs.v;
            T_inlet = inputs.T_inlet;
            T_env = inputs.T_env;
            Q_sinks = inputs.Q_sinks;
        
            T_dot = zeros(N, 1);
        
            % Loop over segments
            for i = 1:N
                T_i = T(i);
                
                % Impose boundary condition on first segment
                if i == 1
                    T_i_minus_1 = T_inlet;
                % Define for other segments
                else
                    T_i_minus_1 = T(i-1);
                end
        
                % Impose boundary condition on last segment
                % i.e zero-gradient at 'outlet'
                if i == N
                    T_i_plus_1 = T_i_minus_1;  
                % Define for other segments
                else
                    T_i_plus_1 = T(i+1);
                end
        
                % Spatial derivatives
                d2T_dx2 = (T_i_plus_1 - 2*T_i + T_i_minus_1) / (Delta_x^2);
        
                if v >= 0
                    dT_dx = (T_i - T_i_minus_1) / Delta_x; % upwind
                else
                    dT_dx = (T_i_plus_1 - T_i) / Delta_x;
                end
        
                % Time derivative
                T_dot(i) = k_prime * d2T_dx2 - v * dT_dx - alpha * (T_i - T_env) ...
                           - Q_sinks(i) / C_seg;
            end
        end

        %% 4. PUMP ELECTRICAL CIRCUIT 
        function Vm_dot = pump_circuit_dynamics(t, Vm, params, inputs)
            
            % Unpacking parameters
            C_motor = params.C_dh_motor;
            R_motor = params.R_dh_motor;
            k_load = params.k_load; 
            
            % Unpacking inputs
            V_in = inputs.V_in;
            Sp_state = inputs.Sp_state; 
            
            if strcmpi(Sp_state, 'CLOSED')
                I_in = (V_in - Vm) / R_motor;
                I_M = k_load * Vm;
                I_C = I_in - I_M;
                Vm_dot = I_C / C_motor;
            else 
                I_M = k_load * Vm;
                I_C = -I_M;
                Vm_dot = I_C / C_motor;
            end
        end
        
        %% 5. HYDROPOWER RESERVOIR 
        function h_dot = reservoir_dynamics(t, h, params, inputs)
            
            if t <= 24*3600
                w_in = inputs.w_in(1);
            else
                w_in = inputs.w_in(2);
            end 

            % Unpacking parameters
            A_res = params.A_res;
            rho_w = params.rho_w;
            g = params.g;
            R_hydro_flow = params.R_hydro_flow;
            h_min = params.h_hydro_min;
            h_max = params.h_hydro_max;
            
            % Unpacking inputs
            % w_in = inputs.w_in;
            
            % Condition on minimum water level
            if h > h_min
                w_out = (rho_w * g * h) / R_hydro_flow;
            else
                w_out = 0;
            end
            
            % Model
            h_dot = (w_in - w_out) / A_res;
            
            if h >= h_max && h_dot > 0
                h_dot = 0; 
            end
        end
        
        %% 6. BACKUP BATTERY 
        function q_dot = battery_dynamics(t, q, params, inputs)
            
            % Unpacking parameters
            R_dh_motor = params.R_dh_motor;
            V_batt = params.V_batt; % 600 V
            I_charging = params.I_charging;
            C_battery = params.C_hydro_battery;
            
            % Unpacking inputs
            Vm = inputs.Vm;
            mode = inputs.mode; % 'DISCHARGE', 'CHARGE', or 'IDLE'

            % Handling switched logic
            if strcmpi(mode, 'DISCHARGE') && q > 0
                I_batt = (V_batt - Vm) / R_dh_motor;
                q_dot = -I_batt;
                
            elseif strcmpi(mode, 'CHARGE') && q < C_battery
                q_dot = I_charging;
                
            else
                % Idle (or full, or empty)
                q_dot = 0; %
            end
        end

    end 
end 
clc, clearvars

%Calculating Reynolds Number at Inlet
velocity_1 = 0.2;
rho = 998.2;
y_plus = 0.95;

pipe_diameter_1 = .06;

Kinematic_viscosity = 1e-6;

Reynolds_Number_1 = (velocity_1*pipe_diameter_1)/Kinematic_viscosity

if Reynolds_Number_1 <= 2300
    flowType_1 = 'Laminar'
elseif Reynolds_Number_1 < 4000 && Reynolds_Number_1 > 2300
       flowType_1 = 'Transition'
else Reynolds_Number_1 >= 4000 
    flowType_1 = 'Turbulent'
end

% Applying Continuity for c = 8

contraction_ratio = 8;
pipe_diameter_2 = pipe_diameter_1 / sqrt(contraction_ratio);
velocity_2 = velocity_1 * contraction_ratio;

%Calculating Reynolds Number at Outlet

Reynolds_Number_2 = (velocity_2*pipe_diameter_2)/Kinematic_viscosity

if Reynolds_Number_2 <= 2300
    flowType_2 = 'Laminar'
elseif Reynolds_Number_2 < 4000 && Reynolds_Number_2 > 2300
       flowType_2 = 'Transition'
else Reynolds_Number_2 >= 4000 
    flowType_2 = 'Turbulent'
end

% Wall y+ Calculator

if Reynolds_Number_1 < 2300
            Cf = 0.664 / sqrt(Re); % Laminar
        else
            Cf = 0.058 * Reynolds_Number_1^(-0.2); % Turbulent
        end
        
        % Wall Shear Stress for Hi
        tau_w_1 = 0.5 * Cf * rho * velocity_1^2;
        
        % Friction Velocity
        u_tau_1= sqrt(tau_w_1 / rho);
        
        % First Layer Height (m)
        Delta_y_1_m = (y_plus * Kinematic_viscosity) / (u_tau_1);
        Delta_y_1_mm = Delta_y_1_m * 1000

        if Reynolds_Number_1 < 2300
            Cf = 0.664 / sqrt(Re); % Laminar
        else
            Cf = 0.058 * Reynolds_Number_1^(-0.2); % Turbulent
        end
    %============================================================        
        
    % Wall Shear Stress for He
        tau_w_2 = 0.5 * Cf * rho * velocity_2^2;
        
        % Friction Velocity
        u_tau_2 = sqrt(tau_w_2 / rho);
        
        % First Layer Height (m)
        Delta_y_2_m = (y_plus * Kinematic_viscosity) / (u_tau_2);
        Delta_y_2_mm = Delta_y_2_m * 1000
clc, clearvars, close all

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

% =========================================================================
%        DIMENSIONLESS VENTURI PLOTTER (Auto-Clean Original CSV)
% =========================================================================

% --- 1. LOAD & CLEAN THE ORIGINAL DATA ---
filename = 'Centreline_velocity_Cleaned_2.csv.txt';

if ~isfile(filename)
    error('File "Values.csv" not found. Make sure it is in the current folder.');
end
% Read the raw text from the file
fid = fopen(filename, 'r');
rawText = fscanf(fid, '%c');
fclose(fid);

% Remove double quotes (") which are causing issues
cleanText = strrep(rawText, '"', '');

% Parse the data (Assume 2 columns: Position, Velocity)
data = sscanf(cleanText, '%f %f'); % Scans all numbers
data = reshape(data, 2, []).';     % Reshape into 2 columns

% Sort data by Position (Column 1) to ensure the line graph is smooth
[~, sortIdx] = sort(data(:, 1));
data = data(sortIdx, :);

% Extract Variables
x_raw = data(:, 1); % Position [m]
V_raw = data(:, 2); % Velocity [m/s]

% --- 2. DEFINE REFERENCE VALUES ---
D_ref = 0.060; 

% Inlet Velocity (V_inlet) - Uses the velocity at the first point
V_ref = V_raw(1); 

% --- 3. CALCULATE DIMENSIONLESS NUMBERS ---
x_star = x_raw / D_ref;      % Dimensionless Position (x / D)
V_star = V_raw / V_ref;      % Velocity Amplification (V / V_inlet)

% --- 4. PLOT WITH WHITE BACKGROUND ---
fig = figure('Color', 'w');  
plot(x_star, V_star, 'LineWidth', 2, 'Color', [0 0.447 0.741]); 

x_throat_actual = (0.6 * D_ref + 0.01) ;  % 
x_throat_dim = x_throat_actual / D_ref; 

% Syntax: xline(value, 'line_style', 'Label', 'Properties'...)
xline(x_throat_dim, ':r', 'Throat', ...
    'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 12, 'FontWeight', 'bold');

% Graph Styling
grid on;
set(gca, 'Color', 'w');             % Set Axis Background to White
set(gca, 'XColor', 'k', 'YColor', 'k'); % Ensure text is Black
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% Labels & Title
xlabel('Non-Dimensional Position (x / D_{inlet})', 'FontWeight', 'bold');
ylabel('Velocity Ratio (V / V_{inlet})', 'FontWeight', 'bold');
title('Dimensionless Centerline Velocity Profile');

% Add reference lines
yline(1, '--k', 'Inlet Velocity', 'LabelHorizontalAlignment', 'left');
xline(0, '--k', 'Inlet Start', 'LabelVerticalAlignment', 'bottom');

axis tight;
clear; clc; close all;

% =========================================================================
% 1. SYSTEM PARAMETERS
% =========================================================================
rho = 998.2;          % Density of Water [kg/m^3]
mu  = 1.003e-3;       % Dynamic Viscosity [Pa.s]
L   = 0.120;          % Diffuser Length [m]
D_inlet = 0.060;      % Inlet Diameter [m]
R_inlet = D_inlet/2;  % Inlet Radius [m]
theta_deg = 5;        % Half Angle [degrees]

% =========================================================================
% 2. READ & PROCESS DATA
% =========================================================================
file_pressure = 'Diffuser_Without_Swirl_CL_Pressure';
file_velocity = 'Diffuser_With_Swirl_CL_Velocity';

[x_p, p_raw] = read_fluent_xy(file_pressure);
[x_v, v_raw] = read_fluent_xy(file_velocity);

% --- Synchronization ---
% Now that x_p and x_v are clean and unique, interp1 will work.
x = x_v; 
v = v_raw;
p = interp1(x_p, p_raw, x, 'linear'); 

% =========================================================================
% 3. CALCULATIONS (Physics & Loss)
% =========================================================================

% A. Geometry Profile (Diameter at every x)
R_local = R_inlet + x .* tand(theta_deg);
D_local = 2 * R_local;

% B. Reynolds Number Profile (Re = rho * v * D / mu)
Re = (rho .* v .* D_local) ./ mu;

% C. Loss Coefficient Calculation
% Calculate Dynamic Pressure (q) and Total Pressure (P_total)
q_local = 0.5 * rho * v.^2;
P_total = p + q_local;

% Define Inlet References (at x = 0)
P_total_in = P_total(1);
q_in       = q_local(1);

% Calculate Cumulative Loss Coefficient
K_loss = (P_total_in - P_total) ./ q_in;

% =========================================================================
% 4. PLOTTING
% =========================================================================

% --- FIGURE 1: Velocity & Pressure (Same Graph) ---
figure('Name', 'Fluid Dynamics', 'Color', 'w', 'Position', [50, 100, 1000, 500]);

yyaxis left
plot(x, p, 'b-', 'LineWidth', 2);
ylabel('Static Pressure [Pa]', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'b'); 

hold on;

yyaxis right
plot(x, v, 'r-', 'LineWidth', 2);
ylabel('Axial Velocity [m/s]', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YColor', 'r');

title('Centerline Distribution: Pressure vs Velocity', 'FontSize', 14);
xlabel('Position [m]', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend({'Static Pressure', 'Axial Velocity'}, 'Location', 'best');


% --- FIGURE 2: Diameter Variation ---
figure('Name', 'Geometry', 'Color', 'w', 'Position', [50, 50, 800, 400]);

plot(x, D_local * 1000, 'k-', 'LineWidth', 2); % Plotting in mm
hold on;
yline(D_inlet*1000, 'k--', 'Inlet Dia');

title('Diameter Variation (Inlet to Outlet)', 'FontSize', 14);
xlabel('Position [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Diameter [mm]', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
area(x, D_local*1000, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none'); % Fill
plot(x, D_local * 1000, 'k-', 'LineWidth', 2); % Redraw line


% --- FIGURE 3: Loss Coefficient vs Reynolds Number ---
figure('Name', 'Loss Analysis', 'Color', 'w', 'Position', [100, 100, 900, 500]);

plot(Re, K_loss, 'm-o', 'LineWidth', 1.5, 'MarkerSize', 3);
grid on;

% Reverse X-Axis (Flow moves from High Re to Low Re)

title(['Loss Coefficient (K) vs Reynolds Number (\theta = ' num2str(theta_deg) '^\circ)'], 'FontSize', 14);
xlabel('Reynolds Number (Re)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Loss Coefficient (K)', 'FontSize', 12, 'FontWeight', 'bold');
subtitle('(X-Axis Reversed: Inlet is on the Left)');

% =========================================================================
% 5. HELPER FUNCTION (UPDATED TO FIX ERRORS)
% =========================================================================
function [x, y] = read_fluent_xy(filename)
    fid = fopen(filename, 'r');
    if fid == -1, error(['Could not open ' filename]); end
    x = []; y = [];
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, '((xy/key/label')
            while ~feof(fid)
                data_line = fgetl(fid);
                % Stop if block ends
                if contains(data_line, ')') && ~any(regexp(data_line, '\d')), break; end
                
                % Clean and parse
                data_line = strrep(data_line, ')', '');
                nums = str2num(data_line);
                if length(nums) == 2, x = [x; nums(1)]; y = [y; nums(2)]; end
            end
            break;
        end
    end
    fclose(fid);
    
    % --- CRITICAL FIX: SORT AND REMOVE DUPLICATES ---
    [x, uniqueIdx] = unique(x); % 'unique' automatically sorts x and removes duplicates
    y = y(uniqueIdx);           % Reorder y to match the sorted x
end
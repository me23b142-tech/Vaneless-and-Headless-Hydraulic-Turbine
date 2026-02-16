clc; clear; close all;

% 1. DEFINE PARAMETERS

Hi = 0.03;              % Inlet Radius [m]
He = 0.03 / sqrt(8);    % Outlet Radius [m] (Area Ratio 8)
L  = 1.2 * Hi;        % Contraction Length [m]

fprintf('--- Design Parameters ---\n');
fprintf('Inlet Radius (Hi): %.4f m\n', Hi);
fprintf('Outlet Radius (He): %.4f m\n', He);
fprintf('Length (L):         %.4f m\n', L);

% 2. SOLVE FOR COEFFICIENTS (Matrix Method)

% Polynomial: Y(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 + c5*x^5
% We solve for [c0, c1, c2, c3, c4, c5] using 6 Boundary Conditions.

% Matrix A (Coefficients of unknown c's)
A = [
    1, 0, 0, 0,      0,      0;       % Y(0)   = Hi
    0, 1, 0, 0,      0,      0;       % Y'(0)  = 0
    0, 0, 2, 0,      0,      0;       % Y''(0) = 0 (Start Curvature)
    1, L, L^2, L^3,   L^4,    L^5;    % Y(L)   = He
    0, 1, 2*L, 3*L^2, 4*L^3,  5*L^4;  % Y'(L)  = 0
    0, 0, 2,   6*L,   12*L^2, 20*L^3  % Y''(L) = 0 (End Curvature)
];

% Vector B (Values)
B = [Hi; 0; 0; He; 0; 0];

% Solve [A]*[C] = [B]
coeffs = A \ B; 

c0 = coeffs(1); c1 = coeffs(2); c2 = coeffs(3);
c3 = coeffs(4); c4 = coeffs(5); c5 = coeffs(6);

% 3. FIND POINT OF INFLECTION (Where Y'' = 0)

% Y''(x) = 2*c2 + 6*c3*x + 12*c4*x^2 + 20*c5*x^3
% We need to find roots of this cubic equation.
% Since c2=0 (from BCs), Y'' = x * (6*c3 + 12*c4*x + 20*c5*x^2)
% One root is x=0. We need the roots of the quadratic part inside brackets.

% Quadratic coefficients for Y'' = A*x^2 + B*x + C
QA = 20 * c5;
QB = 12 * c4;
QC = 6  * c3;

% Roots of quadratic formula
discriminant = sqrt(QB^2 - 4*QA*QC);
root1 = (-QB + discriminant) / (2*QA);
root2 = (-QB - discriminant) / (2*QA);

% Select the valid root within the domain (0 < x < L)
if (root1 > 0 && root1 < L)
    x_inflection = root1;
elseif (root2 > 0 && root2 < L)
    x_inflection = root2;
else
    x_inflection = NaN; % Should not happen for standard Bell-Mehta
end

fprintf('\n--- Inflection Point Analysis ---\n');
fprintf('Calculated Inflection Point x: %.4f m\n', x_inflection);
fprintf('Normalized Location (x/L):     %.2f\n', x_inflection/L);

% 4. PLOTTING

x_plot = linspace(0, L, 500);

% Calculate Y values
Y_plot = c0 + c1*x_plot + c2*x_plot.^2 + c3*x_plot.^3 + c4*x_plot.^4 + c5*x_plot.^5;

% Calculate Y'' values (Curvature proxy)
Y_double_prime = 2*c2 + 6*c3*x_plot + 12*c4*x_plot.^2 + 20*c5*x_plot.^3;

figure('Color','w', 'Position', [100 100 900 600]);

% Subplot 1: The Geometry Profile
subplot(2,1,1);
plot(x_plot, Y_plot, 'b-', 'LineWidth', 2); hold on;
xline(x_inflection, 'k--', 'Inflection');
plot(x_inflection, polyval(flip(coeffs), x_inflection), 'ro', 'MarkerFaceColor','r');
ylabel('Radius Y [m]');
title('Venturi Wall Profile');
grid on;

% Subplot 2: The Curvature (Second Derivative)
subplot(2,1,2);
plot(x_plot, Y_double_prime, 'r-', 'LineWidth', 1.5); hold on;
yline(0, 'k-');
xline(x_inflection, 'k--');
plot(x_inflection, 0, 'bo', 'MarkerFaceColor','b');
xlabel('Axial Distance x [m]');
ylabel('Curvature Proxy (d^2Y/dx^2)');
title('Second Derivative (Showing Inflection Point)');
grid on;

% Ensure data is column-oriented
x_col = 100 * x_plot(:);
y_col = 100 * Y_plot(:);
z_col = zeros(size(x_col)); % Z is 0 for a 2D profile

% --- EXPORT FOR FUSION 360 (.csv) ---
% Format: x, y, z (Unit-less, Fusion usually interprets as cm or mm depending on settings)
% It is safer to convert to the units your CAD is set to. Assuming Meters here:
T_fusion = table(x_col, y_col, z_col);
writetable(T_fusion, 'fusion_curve.csv', 'WriteVariableNames', false);

% --- EXPORT FOR ANSYS DISCOVERY/SPACECLAIM (.txt) ---
% Format: Column 1=Group, Col 2=X, Col 3=Y, Col 4=Z
group_col = ones(size(x_col)); % Group 1 for all points
data_ansys = [group_col, x_col, y_col, z_col];
writematrix(data_ansys, 'ansys_curve.txt', 'Delimiter', 'space');

fclose("all")
fprintf('Files exported: fusion_curve.csv and ansys_curve.txt\n');

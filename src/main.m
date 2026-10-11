% ==============================================================================
% main.m - Main simulation script
% Lactic acid production by Lactiplantibacillus plantarum in a batch
% bioreactor. Tessier kinetic model. Nominal parameters from the assignment.
% Requires: fun_tessier.m
% ==============================================================================

clear all
close all
clc

global Yxs Yps umax Ks kd

% Declare global parameters

%===============================================================================
%                             EXPERIMENTAL DATA
%===============================================================================
% Table 1 - Lactiplantibacillus plantarum producing lactic acid

% Time (h)
texp = [0 1 2 3 4 5 6 7 8 9]';

% Biomass (g/L)
Xexp = [15.0 15.7 17.5 18.6 20.1 21.6 22.4 22.8 22.8 22.9]';

% Substrate - Glucose (g/L)
Sexp = [173.9 159.4 126.8 89.8 57.9 33.6 18.0 6.1 2.5 0.4]';

% Product - Lactic Acid (g/L)
Pexp = [0.0 12.8 26.9 43.3 55.7 67.4 73.2 78.7 81.4 82.7]';

%===============================================================================
%                                TIME (HOURS)
%===============================================================================
ti = 0;
dt = 1;        % Time step in hours
tf = 9;
t = [ti:dt:tf]';

%===============================================================================
%                            INITIAL CONDITIONS
%===============================================================================
y0 = [15.0 173.9 0.0];  % Initial condition [X0 S0 P0]

%===============================================================================
%                             PARAMETER VALUES
%===============================================================================
Yxs = 0.04;         % Biomass yield on substrate (g biomass/g substrate)
Yps = 0.43;         % Product yield on substrate (g product/g substrate)
umax = 0.57;        % Maximum specific growth rate (h^-1)
Ks = 21;            % Saturation constant (g/L)
kd = 0.3;           % Cell death rate (h^-1)

%===============================================================================
%                         SYSTEM OF ODEs SOLUTION
%===============================================================================
[tsim, y] = ode15s(@fun_tessier, t, y0);

Xpred = y(:, 1);    % Extract biomass predictions
Spred = y(:, 2);    % Extract substrate predictions
Ppred = y(:, 3);    % Extract product predictions

%===============================================================================
%                      ERROR CALCULATION AND METRICS
%===============================================================================
% Linear interpolation of experimental data to simulation time points
Xexp_interp = interp1(texp, Xexp, tsim, 'linear');
Sexp_interp = interp1(texp, Sexp, tsim, 'linear');
Pexp_interp = interp1(texp, Pexp, tsim, 'linear');

% Root Mean Square Error (RMSE)
RMSE_X = sqrt(mean((Xpred - Xexp_interp).^2));
RMSE_S = sqrt(mean((Spred - Sexp_interp).^2));
RMSE_P = sqrt(mean((Ppred - Pexp_interp).^2));

fprintf('\n========== ERROR METRICS ==========\n');
fprintf('RMSE Biomass (X): %.4f g/L\n', RMSE_X);
fprintf('RMSE Substrate (S): %.4f g/L\n', RMSE_S);
fprintf('RMSE Product (P): %.4f g/L\n', RMSE_P);
fprintf('===================================\n\n');

%===============================================================================
%                                PLOTS
%===============================================================================

% Figure 1: Biomass
figure(1)
plot(tsim, Xpred, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Predicted');
hold on
plot(texp, Xexp, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'Experimental');
xlabel('\bf \fontname{Times New Roman} Time (h)', 'Fontsize', 28);
ylabel('\bf \fontname{Times New Roman} X (g/L)', 'Fontsize', 28);
title('\bf \fontname{Times New Roman} Biomass Growth', 'Fontsize', 32);
legend('Location', 'best', 'Fontsize', 24);
grid on
set(gca, 'Fontsize', 24, 'fontWeight', 'bold', 'fontname', 'Times New Roman');
set(gcf, 'Position', [100 400 900 600]);

% Figure 2: Substrate
figure(2)
plot(tsim, Spred, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Predicted');
hold on
plot(texp, Sexp, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Experimental');
xlabel('\bf \fontname{Times New Roman} Time (h)', 'Fontsize', 28);
ylabel('\bf \fontname{Times New Roman} S (g/L)', 'Fontsize', 28);
title('\bf \fontname{Times New Roman} Substrate Consumption (Glucose)', 'Fontsize', 32);
legend('Location', 'best', 'Fontsize', 24);
grid on
set(gca, 'Fontsize', 24, 'fontWeight', 'bold', 'fontname', 'Times New Roman');
set(gcf, 'Position', [1050 400 900 600]);

% Figure 3: Product
figure(3)
plot(tsim, Ppred, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Predicted');
hold on
plot(texp, Pexp, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 10, 'DisplayName', 'Experimental');
xlabel('\bf \fontname{Times New Roman} Time (h)', 'Fontsize', 28);
ylabel('\bf \fontname{Times New Roman} P (g/L)', 'Fontsize', 28);
title('\bf \fontname{Times New Roman} Lactic Acid Production', 'Fontsize', 32);
legend('Location', 'best', 'Fontsize', 24);
grid on
set(gca, 'Fontsize', 24, 'fontWeight', 'bold', 'fontname', 'Times New Roman');
set(gcf, 'Position', [2000 400 900 600]);

% Figure 4: Combined Plot
figure(4)
subplot(1, 3, 1)
plot(tsim, Xpred, 'b-', 'LineWidth', 2.5);
hold on
plot(texp, Xexp, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('\bf Time (h)', 'Fontsize', 20);
ylabel('\bf X (g/L)', 'Fontsize', 20);
title('\bf Biomass', 'Fontsize', 22, 'fontweight', 'bold');
grid on
legend('Predicted', 'Experimental', 'Fontsize', 16);
set(gca, 'Fontsize', 18);

subplot(1, 3, 2)
plot(tsim, Spred, 'b-', 'LineWidth', 2.5);
hold on
plot(texp, Sexp, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
xlabel('\bf Time (h)', 'Fontsize', 20);
ylabel('\bf S (g/L)', 'Fontsize', 20);
title('\bf Substrate', 'Fontsize', 22, 'fontweight', 'bold');
grid on
legend('Predicted', 'Experimental', 'Fontsize', 16);
set(gca, 'Fontsize', 18);

subplot(1, 3, 3)
plot(tsim, Ppred, 'b-', 'LineWidth', 2.5);
hold on
plot(texp, Pexp, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
xlabel('\bf Time (h)', 'Fontsize', 20);
ylabel('\bf P (g/L)', 'Fontsize', 20);
title('\bf Product', 'Fontsize', 22, 'fontweight', 'bold');
grid on
legend('Predicted', 'Experimental', 'Fontsize', 16);
set(gca, 'Fontsize', 18);

set(gcf, 'Position', [100 50 1400 600]);

fprintf('\nSimulation successfully completed!\n');
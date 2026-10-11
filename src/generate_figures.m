clear all
close all
clc

global Yxs Yps umax Ks kd

% ==============================================================================
% Generates the article figures in PDF format (folder ../figures)
% ==============================================================================

folder = '../figures';
if ~exist(folder, 'dir')
    mkdir(folder);
end

%===============================================================================
% Experimental data and parameters
%===============================================================================
texp = [0 1 2 3 4 5 6 7 8 9]';
Xexp = [15.0 15.7 17.5 18.6 20.1 21.6 22.4 22.8 22.8 22.9]';
Sexp = [173.9 159.4 126.8 89.8 57.9 33.6 18.0 6.1 2.5 0.4]';
Pexp = [0.0 12.8 26.9 43.3 55.7 67.4 73.2 78.7 81.4 82.7]';
y0 = [15.0 173.9 0.0];

tfine = linspace(0, 9, 200)';   % fine grid for smooth curves

Yxs = 0.04; Yps = 0.43; Ks = 21;

% Simulation with nominal parameters
umax = 0.57; kd = 0.30;
[~, y_nom] = ode15s(@fun_tessier, tfine, y0);

% Simulation with adjusted parameters
umax = 0.0705; kd = 0.0;
[~, y_adj] = ode15s(@fun_tessier, tfine, y0);
[~, y_adj_p] = ode15s(@fun_tessier, texp, y0);   % at experimental sampling points

%===============================================================================
% Figure 1 - Fit of the three state variables
%===============================================================================
f1 = figure('Position', [100 100 1500 450]);

data = [Xexp Sexp Pexp];
lab  = {'X (g L^{-1})', 'S (g L^{-1})', 'P (g L^{-1})'};
tit  = {'(a) Biomass', '(b) Substrate', '(c) Lactic acid'};

for k = 1:3
    subplot(1, 3, k)
    plot(tfine, y_adj(:,k), 'k-', 'LineWidth', 2.2)
    hold on
    plot(tfine, y_nom(:,k), 'k--', 'LineWidth', 1.6, 'Color', [0.55 0.55 0.55])
    plot(texp, data(:,k), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.4)
    xlabel('Time (h)', 'FontSize', 13)
    ylabel(lab{k}, 'FontSize', 13)
    title(tit{k}, 'FontSize', 14)
    xlim([0 9])
    grid on
    set(gca, 'FontSize', 12, 'GridAlpha', 0.25)
    if k == 1
        legend({'Fitted model', 'Nominal parameters', 'Experimental'}, ...
               'Location', 'southeast', 'FontSize', 10);
    end
end

print(f1, fullfile(folder, 'fig_fit.pdf'), '-dpdf', '-S1500,450');

%===============================================================================
% Figure 2 - Specific growth rate as a function of substrate concentration
%===============================================================================
f2 = figure('Position', [100 100 750 500]);

Sv = linspace(0, 180, 400);
mu_tes_adj  = 0.0705 * (1 - exp(-Sv / 21));
mu_tes_nom  = 0.5700 * (1 - exp(-Sv / 21));
mu_monod    = 0.0705 * Sv ./ (21 + Sv);

plot(Sv, mu_tes_nom, '--', 'LineWidth', 1.8, 'Color', [0.55 0.55 0.55])
hold on
plot(Sv, mu_tes_adj, 'k-', 'LineWidth', 2.2)
plot(Sv, mu_monod, 'k-.', 'LineWidth', 1.6)

% Substrate range covered by the experiment
yl = ylim;
plot([min(Sexp) min(Sexp)], [0 yl(2)], ':', 'Color', [0.3 0.3 0.3])
plot([max(Sexp) max(Sexp)], [0 yl(2)], ':', 'Color', [0.3 0.3 0.3])

xlabel('S (g L^{-1})', 'FontSize', 13)
ylabel('\mu (h^{-1})', 'FontSize', 13)
legend({'Tessier (nominal \mu_{max} = 0.57 h^{-1})', ...
        'Tessier (fitted \mu_{max} = 0.0705 h^{-1})', ...
        'Monod (fitted \mu_{max}, same K_S)'}, ...
       'Location', 'east', 'FontSize', 10)
grid on
set(gca, 'FontSize', 12, 'GridAlpha', 0.25)

print(f2, fullfile(folder, 'fig_mu.pdf'), '-dpdf', '-S750,500');

%===============================================================================
% Figure 3 - Parity plot (predicted vs. experimental)
%===============================================================================
f3 = figure('Position', [100 100 600 560]);

plot(Xexp, y_adj_p(:,1), 'ks', 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 8)
hold on
plot(Sexp, y_adj_p(:,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.3)
plot(Pexp, y_adj_p(:,3), 'k^', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 8)
plot([0 180], [0 180], 'k-', 'LineWidth', 1.2)
plot([0 180], [0 180]*1.15, 'k:', 'LineWidth', 1)
plot([0 180], [0 180]*0.85, 'k:', 'LineWidth', 1)

xlabel('Experimental value (g L^{-1})', 'FontSize', 13)
ylabel('Predicted value (g L^{-1})', 'FontSize', 13)
legend({'X', 'S', 'P', 'y = x', '\pm 15 %'}, 'Location', 'southeast', 'FontSize', 11)
axis([0 180 0 180])
axis square
grid on
set(gca, 'FontSize', 12, 'GridAlpha', 0.25)

print(f3, fullfile(folder, 'fig_parity.pdf'), '-dpdf', '-S600,560');

printf('Figures generated in %s\n', folder);
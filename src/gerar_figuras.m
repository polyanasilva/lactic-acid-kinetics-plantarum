clear all
close all
clc

global Yxs Yps umax Ks kd

% ==============================================================================
% Gera as figuras do artigo em formato PDF (pasta ../figuras)
% ==============================================================================

pasta = '../figuras';
if ~exist(pasta, 'dir')
    mkdir(pasta);
end

%===============================================================================
% Dados experimentais e parâmetros
%===============================================================================
texp = [0 1 2 3 4 5 6 7 8 9]';
Xexp = [15.0 15.7 17.5 18.6 20.1 21.6 22.4 22.8 22.8 22.9]';
Sexp = [173.9 159.4 126.8 89.8 57.9 33.6 18.0 6.1 2.5 0.4]';
Pexp = [0.0 12.8 26.9 43.3 55.7 67.4 73.2 78.7 81.4 82.7]';
y0 = [15.0 173.9 0.0];

tfino = linspace(0, 9, 200)';   % malha fina para curvas suaves

Yxs = 0.04; Yps = 0.43; Ks = 21;

% Simulação com parâmetros nominais
umax = 0.57; kd = 0.30;
[~, y_nom] = ode15s(@fun_tessier, tfino, y0);

% Simulação com parâmetros ajustados
umax = 0.0705; kd = 0.0;
[~, y_aj] = ode15s(@fun_tessier, tfino, y0);
[~, y_aj_p] = ode15s(@fun_tessier, texp, y0);   % nos pontos experimentais

%===============================================================================
% Figura 1 - Ajuste das três variáveis de estado
%===============================================================================
f1 = figure('Position', [100 100 1500 450]);

dados = [Xexp Sexp Pexp];
lab   = {'X (g L^{-1})', 'S (g L^{-1})', 'P (g L^{-1})'};
tit   = {'(a) Biomassa', '(b) Substrato', '(c) Ácido lático'};

for k = 1:3
    subplot(1, 3, k)
    plot(tfino, y_aj(:,k), 'k-', 'LineWidth', 2.2)
    hold on
    plot(tfino, y_nom(:,k), 'k--', 'LineWidth', 1.6, 'Color', [0.55 0.55 0.55])
    plot(texp, dados(:,k), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.4)
    xlabel('Tempo (h)', 'FontSize', 13)
    ylabel(lab{k}, 'FontSize', 13)
    title(tit{k}, 'FontSize', 14)
    xlim([0 9])
    grid on
    set(gca, 'FontSize', 12, 'GridAlpha', 0.25)
    if k == 1
        legend({'Modelo ajustado', 'Parâmetros nominais', 'Experimental'}, ...
               'Location', 'southeast', 'FontSize', 10);
    end
end

print(f1, fullfile(pasta, 'fig_ajuste.pdf'), '-dpdf', '-S1500,450');

%===============================================================================
% Figura 2 - Taxa específica de crescimento em função do substrato
%===============================================================================
f2 = figure('Position', [100 100 750 500]);

Sv = linspace(0, 180, 400);
mu_tes_aj  = 0.0705 * (1 - exp(-Sv / 21));
mu_tes_nom = 0.5700 * (1 - exp(-Sv / 21));
mu_monod   = 0.0705 * Sv ./ (21 + Sv);

plot(Sv, mu_tes_nom, '--', 'LineWidth', 1.8, 'Color', [0.55 0.55 0.55])
hold on
plot(Sv, mu_tes_aj, 'k-', 'LineWidth', 2.2)
plot(Sv, mu_monod, 'k-.', 'LineWidth', 1.6)

% Faixa de substrato coberta pelo experimento
yl = ylim;
plot([min(Sexp) min(Sexp)], [0 yl(2)], ':', 'Color', [0.3 0.3 0.3])
plot([max(Sexp) max(Sexp)], [0 yl(2)], ':', 'Color', [0.3 0.3 0.3])

xlabel('S (g L^{-1})', 'FontSize', 13)
ylabel('\mu (h^{-1})', 'FontSize', 13)
legend({'Tessier (\mu_{max} nominal = 0,57 h^{-1})', ...
        'Tessier (\mu_{max} ajustado = 0,0705 h^{-1})', ...
        'Monod (\mu_{max} ajustado, mesmo K_S)'}, ...
       'Location', 'east', 'FontSize', 10)
grid on
set(gca, 'FontSize', 12, 'GridAlpha', 0.25)

print(f2, fullfile(pasta, 'fig_mu.pdf'), '-dpdf', '-S750,500');

%===============================================================================
% Figura 3 - Gráfico de paridade (predito vs. experimental)
%===============================================================================
f3 = figure('Position', [100 100 600 560]);

plot(Xexp, y_aj_p(:,1), 'ks', 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 8)
hold on
plot(Sexp, y_aj_p(:,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.3)
plot(Pexp, y_aj_p(:,3), 'k^', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 8)
plot([0 180], [0 180], 'k-', 'LineWidth', 1.2)
plot([0 180], [0 180]*1.15, 'k:', 'LineWidth', 1)
plot([0 180], [0 180]*0.85, 'k:', 'LineWidth', 1)

xlabel('Valor experimental (g L^{-1})', 'FontSize', 13)
ylabel('Valor predito (g L^{-1})', 'FontSize', 13)
legend({'X', 'S', 'P', 'y = x', '\pm 15 %'}, 'Location', 'southeast', 'FontSize', 11)
axis([0 180 0 180])
axis square
grid on
set(gca, 'FontSize', 12, 'GridAlpha', 0.25)

print(f3, fullfile(pasta, 'fig_paridade.pdf'), '-dpdf', '-S600,560');

printf('Figuras geradas em %s\n', pasta);

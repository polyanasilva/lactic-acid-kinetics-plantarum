% ==============================================================================
% main.m - Roteiro principal de simulacao
% Producao de acido latico por Lactiplantibacillus plantarum em biorreator
% batelada. Modelo cinetico de Tessier. Parametros nominais do enunciado.
% Requer: fun_tessier.m
% ==============================================================================

clear all
close all
clc

global Yxs Yps umax Ks kd

% Declara parâmetros globais

%===============================================================================
%                             DADOS EXPERIMENTAIS
%===============================================================================
% Tabela 1 - Lactiplantibacillus plantarum produzindo ácido lático

% Tempo (h)
texp = [0 1 2 3 4 5 6 7 8 9]';

% Biomassa (g/L)
Xexp = [15.0 15.7 17.5 18.6 20.1 21.6 22.4 22.8 22.8 22.9]';

% Substrato - Glicose (g/L)
Sexp = [173.9 159.4 126.8 89.8 57.9 33.6 18.0 6.1 2.5 0.4]';

% Produto - Ácido Lático (g/L)
Pexp = [0.0 12.8 26.9 43.3 55.7 67.4 73.2 78.7 81.4 82.7]';

%===============================================================================
%                             TEMPO (HORAS)
%===============================================================================
ti = 0;
dt = 1;        % Intervalo de tempo em horas
tf = 9;
t = [ti:dt:tf]';

%===============================================================================
%                            CONDIÇÃO INICIAL
%===============================================================================
y0 = [15.0 173.9 0.0];  % Condição inicial [X0 S0 P0]

%===============================================================================
%                            VALOR DOS PARÂMETROS
%===============================================================================
Yxs = 0.04;         % Rendimento de biomassa em relação ao substrato (g biomassa/g substrato)
Yps = 0.43;         % Rendimento de produto em relação ao substrato (g produto/g substrato)
umax = 0.57;        % Taxa máxima de crescimento específico (h^-1)
Ks = 21;            % Constante de saturação (g/L)
kd = 0.3;           % Taxa de morte celular (h^-1)

%===============================================================================
%                       SOLUÇÃO DO SISTEMA DE EDOs
%===============================================================================
[tsim, y] = ode15s(@fun_tessier, t, y0);

Xpred = y(:, 1);    % Extrai as predições de biomassa
Spred = y(:, 2);    % Extrai as predições de substrato
Ppred = y(:, 3);    % Extrai as predições de produto

%===============================================================================
%                      CÁLCULO DE ERROS E MÉTRICAS
%===============================================================================
% Interpolação dos dados experimentais para os tempos simulados
Xexp_interp = interp1(texp, Xexp, tsim, 'linear');
Sexp_interp = interp1(texp, Sexp, tsim, 'linear');
Pexp_interp = interp1(texp, Pexp, tsim, 'linear');

% Erro médio quadrático (RMSE)
RMSE_X = sqrt(mean((Xpred - Xexp_interp).^2));
RMSE_S = sqrt(mean((Spred - Sexp_interp).^2));
RMSE_P = sqrt(mean((Ppred - Pexp_interp).^2));

fprintf('\n========== MÉTRICAS DE ERRO ==========\n');
fprintf('RMSE Biomassa (X): %.4f g/L\n', RMSE_X);
fprintf('RMSE Substrato (S): %.4f g/L\n', RMSE_S);
fprintf('RMSE Produto (P): %.4f g/L\n', RMSE_P);
fprintf('======================================\n\n');

%===============================================================================
%                              GRÁFICOS
%===============================================================================

% Figura 1: Biomassa
figure(1)
plot(tsim, Xpred, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Predito');
hold on
plot(texp, Xexp, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'Experimental');
xlabel('\bf \fontname{Times New Roman} Tempo (h)', 'Fontsize', 28);
ylabel('\bf \fontname{Times New Roman} X (g/L)', 'Fontsize', 28);
title('\bf \fontname{Times New Roman} Crescimento de Biomassa', 'Fontsize', 32);
legend('Location', 'best', 'Fontsize', 24);
grid on
set(gca, 'Fontsize', 24, 'fontWeight', 'bold', 'fontname', 'Times New Roman');
set(gcf, 'Position', [100 400 900 600]);

% Figura 2: Substrato
figure(2)
plot(tsim, Spred, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Predito');
hold on
plot(texp, Sexp, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Experimental');
xlabel('\bf \fontname{Times New Roman} Tempo (h)', 'Fontsize', 28);
ylabel('\bf \fontname{Times New Roman} S (g/L)', 'Fontsize', 28);
title('\bf \fontname{Times New Roman} Consumo de Substrato (Glicose)', 'Fontsize', 32);
legend('Location', 'best', 'Fontsize', 24);
grid on
set(gca, 'Fontsize', 24, 'fontWeight', 'bold', 'fontname', 'Times New Roman');
set(gcf, 'Position', [1050 400 900 600]);

% Figura 3: Produto
figure(3)
plot(tsim, Ppred, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Predito');
hold on
plot(texp, Pexp, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 10, 'DisplayName', 'Experimental');
xlabel('\bf \fontname{Times New Roman} Tempo (h)', 'Fontsize', 28);
ylabel('\bf \fontname{Times New Roman} P (g/L)', 'Fontsize', 28);
title('\bf \fontname{Times New Roman} Produção de Ácido Lático', 'Fontsize', 32);
legend('Location', 'best', 'Fontsize', 24);
grid on
set(gca, 'Fontsize', 24, 'fontWeight', 'bold', 'fontname', 'Times New Roman');
set(gcf, 'Position', [2000 400 900 600]);

% Figura 4: Gráfico combinado
figure(4)
subplot(1, 3, 1)
plot(tsim, Xpred, 'b-', 'LineWidth', 2.5);
hold on
plot(texp, Xexp, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('\bf Tempo (h)', 'Fontsize', 20);
ylabel('\bf X (g/L)', 'Fontsize', 20);
title('\bf Biomassa', 'Fontsize', 22, 'fontweight', 'bold');
grid on
legend('Predito', 'Experimental', 'Fontsize', 16);
set(gca, 'Fontsize', 18);

subplot(1, 3, 2)
plot(tsim, Spred, 'b-', 'LineWidth', 2.5);
hold on
plot(texp, Sexp, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
xlabel('\bf Tempo (h)', 'Fontsize', 20);
ylabel('\bf S (g/L)', 'Fontsize', 20);
title('\bf Substrato', 'Fontsize', 22, 'fontweight', 'bold');
grid on
legend('Predito', 'Experimental', 'Fontsize', 16);
set(gca, 'Fontsize', 18);

subplot(1, 3, 3)
plot(tsim, Ppred, 'b-', 'LineWidth', 2.5);
hold on
plot(texp, Pexp, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
xlabel('\bf Tempo (h)', 'Fontsize', 20);
ylabel('\bf P (g/L)', 'Fontsize', 20);
title('\bf Produto', 'Fontsize', 22, 'fontweight', 'bold');
grid on
legend('Predito', 'Experimental', 'Fontsize', 16);
set(gca, 'Fontsize', 18);

set(gcf, 'Position', [100 50 1400 600]);

fprintf('\nSimulação concluída com sucesso!\n');

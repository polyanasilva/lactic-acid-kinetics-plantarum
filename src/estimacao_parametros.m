clear all
close all
clc

global Yxs Yps umax Ks kd

% ==============================================================================
% ESTIMAÇÃO DE PARÂMETROS - MODELO DE TESSIER (equação 14)
% Lactiplantibacillus plantarum - produção de ácido lático
%
% Objetivo: reestimar os parâmetros cinéticos por ajuste aos dados
% experimentais da Tabela 1, partindo dos valores nominais do enunciado.
% ==============================================================================

%===============================================================================
%    FUNÇÕES AUXILIARES (no Octave, devem ser definidas antes do uso)
%===============================================================================
% Função objetivo: soma dos quadrados dos resíduos normalizados
function e = funcao_objetivo(p, texp, Xexp, Sexp, Pexp, y0)
    global umax kd
    p = abs(p);
    umax = p(1);
    kd   = p(2);
    try
        [~, y] = ode15s(@fun_tessier, texp, y0);
        if numel(y(:,1)) != numel(texp)
            e = 1e8;   % integração falhou/incompleta: penaliza
            return
        end
        e = sum(((y(:,1) - Xexp) / max(Xexp)).^2) + ...
            sum(((y(:,2) - Sexp) / max(Sexp)).^2) + ...
            sum(((y(:,3) - Pexp) / max(Pexp)).^2);
    catch
        e = 1e8;
    end
end

% Métricas de qualidade do ajuste: RMSE por variável e R2 global
function [rX, rS, rP, R2] = metricas(y, Xexp, Sexp, Pexp)
    rX = sqrt(mean((y(:,1) - Xexp).^2));
    rS = sqrt(mean((y(:,2) - Sexp).^2));
    rP = sqrt(mean((y(:,3) - Pexp).^2));

    SSE = sum((y(:,1) - Xexp).^2) + sum((y(:,2) - Sexp).^2) + sum((y(:,3) - Pexp).^2);
    SST = sum((Xexp - mean(Xexp)).^2) + sum((Sexp - mean(Sexp)).^2) + ...
          sum((Pexp - mean(Pexp)).^2);
    R2 = 1 - SSE / SST;
end

%===============================================================================
%                             DADOS EXPERIMENTAIS
%===============================================================================
texp = [0 1 2 3 4 5 6 7 8 9]';
Xexp = [15.0 15.7 17.5 18.6 20.1 21.6 22.4 22.8 22.8 22.9]';
Sexp = [173.9 159.4 126.8 89.8 57.9 33.6 18.0 6.1 2.5 0.4]';
Pexp = [0.0 12.8 26.9 43.3 55.7 67.4 73.2 78.7 81.4 82.7]';

y0 = [15.0 173.9 0.0];   % Condição inicial [X0 S0 P0]

%===============================================================================
%              PARÂMETROS FIXOS (mantidos do enunciado)
%===============================================================================
% Yxs e Yps são coeficientes estequiométricos de rendimento, coerentes com o
% balanço global dos dados: Yxs = DX/DS = 7,9/173,5 = 0,0455
%                           Yps = DP/DS = 82,7/173,5 = 0,477
Yxs = 0.04;    % g biomassa / g substrato
Yps = 0.43;    % g produto  / g substrato

% Ks é mantido em 21 g/L (valor do enunciado). A análise de sensibilidade ao
% final deste script mostra que Ks é pouco identificável com estes dados,
% pois S >> Ks durante quase toda a batelada (Tessier opera saturado).
Ks = 21;       % g/L

%===============================================================================
%              SIMULAÇÃO COM OS PARÂMETROS NOMINAIS DO ENUNCIADO
%===============================================================================
umax = 0.57;   % h^-1
kd   = 0.3;    % h^-1

[t_nom, y_nom] = ode15s(@fun_tessier, texp, y0);
[rX_nom, rS_nom, rP_nom, R2_nom] = metricas(y_nom, Xexp, Sexp, Pexp);

%===============================================================================
%                        ESTIMAÇÃO DE PARÂMETROS
%===============================================================================
% Parâmetros estimados: umax e kd
% Método: minimização da soma dos quadrados dos resíduos normalizados
%         (Nelder-Mead simplex, função fminsearch do Octave)
%
% A normalização por max(.) de cada variável equaliza o peso de X, S e P na
% função objetivo, evitando que S (ordem de 174 g/L) domine sobre X (ordem
% de 23 g/L).

p0 = [0.10 0.05];   % Estimativa inicial [umax kd]

opcoes = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, ...
                  'TolX', 1e-10, 'TolFun', 1e-10);

p_est = fminsearch(@(p) funcao_objetivo(p, texp, Xexp, Sexp, Pexp, y0), p0, opcoes);
p_est = abs(p_est);   % garante positividade dos parâmetros

umax = p_est(1);
kd   = p_est(2);

[t_aj, y_aj] = ode15s(@fun_tessier, texp, y0);
[rX_aj, rS_aj, rP_aj, R2_aj] = metricas(y_aj, Xexp, Sexp, Pexp);

%===============================================================================
%                          RELATÓRIO DE RESULTADOS
%===============================================================================
printf('\n');
printf('===============================================================\n');
printf('        ESTIMACAO DE PARAMETROS - MODELO DE TESSIER\n');
printf('===============================================================\n');
printf('Parametro          Enunciado        Ajustado\n');
printf('---------------------------------------------------------------\n');
printf('umax (1/h)         %9.4f      %9.4f\n', 0.57, umax);
printf('Ks   (g/L)         %9.4f      %9.4f   (fixo)\n', 21, Ks);
printf('kd   (1/h)         %9.4f      %9.4f\n', 0.30, kd);
printf('Yxs  (g/g)         %9.4f      %9.4f   (fixo)\n', Yxs, Yxs);
printf('Yps  (g/g)         %9.4f      %9.4f   (fixo)\n', Yps, Yps);
printf('---------------------------------------------------------------\n');
printf('RMSE X (g/L)       %9.2f      %9.2f\n', rX_nom, rX_aj);
printf('RMSE S (g/L)       %9.2f      %9.2f\n', rS_nom, rS_aj);
printf('RMSE P (g/L)       %9.2f      %9.2f\n', rP_nom, rP_aj);
printf('R2 global          %9.4f      %9.4f\n', R2_nom, R2_aj);
printf('===============================================================\n\n');

printf('Comparacao ponto a ponto (modelo ajustado):\n');
printf('  t(h)    Xsim    Xexp     Ssim    Sexp     Psim    Pexp\n');
printf('%6.1f %7.2f %7.2f %8.2f %7.2f %8.2f %7.2f\n', ...
       [texp y_aj(:,1) Xexp y_aj(:,2) Sexp y_aj(:,3) Pexp]');
printf('\n');

%===============================================================================
%                   ANÁLISE DE SENSIBILIDADE AO PARÂMETRO Ks
%===============================================================================
% Justifica a decisão de manter Ks = 21 g/L: para cada Ks fixado, reestima-se
% umax e verifica-se o R2 resultante.

printf('Analise de sensibilidade a Ks (umax reestimado para cada Ks):\n');
printf('   Ks (g/L)     umax (1/h)      R2\n');
Ks_teste = [5 10 21 30 50 100];
for i = 1:length(Ks_teste)
    Ks = Ks_teste(i);
    kd = 0;
    p = abs(fminsearch(@(p) funcao_objetivo([p 0], texp, Xexp, Sexp, Pexp, y0), 0.07, opcoes));
    umax = p;
    [~, y_s] = ode15s(@fun_tessier, texp, y0);
    [~, ~, ~, R2_s] = metricas(y_s, Xexp, Sexp, Pexp);
    printf('%10.1f %14.5f %9.5f\n', Ks, umax, R2_s);
end
printf('\n');

% Restaura os parâmetros ajustados para os gráficos
Ks = 21; umax = p_est(1); kd = p_est(2);

%===============================================================================
%                                GRÁFICOS
%===============================================================================
tit = {'Biomassa - X', 'Substrato (glicose) - S', 'Ácido lático - P'};
lab = {'X (g/L)', 'S (g/L)', 'P (g/L)'};
dad = [Xexp Sexp Pexp];

for k = 1:3
    figure(k)
    plot(t_aj, y_aj(:,k), 'b-', 'LineWidth', 2.5)
    hold on
    plot(t_nom, y_nom(:,k), 'r--', 'LineWidth', 2)
    plot(texp, dad(:,k), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 9)
    xlabel('\bf \fontname{Times New Roman} Tempo (h)', 'Fontsize', 24)
    ylabel(['\bf \fontname{Times New Roman} ' lab{k}], 'Fontsize', 24)
    title(['\bf \fontname{Times New Roman} ' tit{k}], 'Fontsize', 26)
    legend({'Ajustado', 'Enunciado', 'Experimental'}, 'Fontsize', 18)
    grid on
    set(gca, 'Fontsize', 18, 'fontWeight', 'bold', 'fontname', 'Times New Roman')
end

printf('Estimacao concluida.\n');

clear all
close all
clc

global Yxs Yps umax Ks kd

% ==============================================================================
% PARAMETER ESTIMATION - TESSIER MODEL (equation 14)
% Lactiplantibacillus plantarum - lactic acid production
%
% Objective: Reestimate the kinetic parameters by fitting to the experimental
% data from Table 1, starting from the nominal values given in the assignment.
% ==============================================================================

%===============================================================================
%    AUXILIARY FUNCTIONS (in Octave, must be defined before use)
%===============================================================================
% Objective function: sum of squared normalized residuals
function e = funcao_objetivo(p, texp, Xexp, Sexp, Pexp, y0)
    global umax kd
    p = abs(p);
    umax = p(1);
    kd   = p(2);
    try
        [~, y] = ode15s(@fun_tessier, texp, y0);
        if numel(y(:,1)) != numel(texp)
            e = 1e8;   % integration failed/incomplete: penalize
            return
        end
        e = sum(((y(:,1) - Xexp) / max(Xexp)).^2) + ...
            sum(((y(:,2) - Sexp) / max(Sexp)).^2) + ...
            sum(((y(:,3) - Pexp) / max(Pexp)).^2);
    catch
        e = 1e8;
    end
end

% Fit quality metrics: RMSE per variable and global R2
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
%                             EXPERIMENTAL DATA
%===============================================================================
texp = [0 1 2 3 4 5 6 7 8 9]';
Xexp = [15.0 15.7 17.5 18.6 20.1 21.6 22.4 22.8 22.8 22.9]';
Sexp = [173.9 159.4 126.8 89.8 57.9 33.6 18.0 6.1 2.5 0.4]';
Pexp = [0.0 12.8 26.9 43.3 55.7 67.4 73.2 78.7 81.4 82.7]';

y0 = [15.0 173.9 0.0];   % Initial condition [X0 S0 P0]

%===============================================================================
%              FIXED PARAMETERS (kept from the assignment)
%===============================================================================
% Yxs and Yps are stoichiometric yield coefficients, consistent with the
% global data balance: Yxs = DX/DS = 7.9/173.5 = 0.0455
%                      Yps = DP/DS = 82.7/173.5 = 0.477
Yxs = 0.04;    % g biomass / g substrate
Yps = 0.43;    % g product / g substrate

% Ks is kept at 21 g/L (assignment value). The sensitivity analysis at the
% end of this script shows that Ks is poorly identifiable with these data,
% as S >> Ks throughout almost the entire batch (Tessier operates saturated).
Ks = 21;       % g/L

%===============================================================================
%              SIMULATION WITH NOMINAL PARAMETERS FROM ASSIGNMENT
%===============================================================================
umax = 0.57;   % h^-1
kd   = 0.3;    % h^-1

[t_nom, y_nom] = ode15s(@fun_tessier, texp, y0);
[rX_nom, rS_nom, rP_nom, R2_nom] = metricas(y_nom, Xexp, Sexp, Pexp);

%===============================================================================
%                        PARAMETER ESTIMATION
%===============================================================================
% Estimated parameters: umax and kd
% Method: Minimization of the sum of squared normalized residuals
%         (Nelder-Mead simplex, Octave's fminsearch function)
%
% Normalization by max(.) of each variable equalizes the weight of X, S, and P
% in the objective function, preventing S (order of 174 g/L) from dominating
% over X (order of 23 g/L).

p0 = [0.10 0.05];   % Initial guess [umax kd]

options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, ...
                   'TolX', 1e-10, 'TolFun', 1e-10);

p_est = fminsearch(@(p) funcao_objetivo(p, texp, Xexp, Sexp, Pexp, y0), p0, options);
p_est = abs(p_est);   % ensures parameter positivity

umax = p_est(1);
kd   = p_est(2);

[t_aj, y_aj] = ode15s(@fun_tessier, texp, y0);
[rX_aj, rS_aj, rP_aj, R2_aj] = metricas(y_aj, Xexp, Sexp, Pexp);

%===============================================================================
%                          RESULTS REPORT
%===============================================================================
printf('\n');
printf('===============================================================\n');
printf('        PARAMETER ESTIMATION - TESSIER MODEL\n');
printf('===============================================================\n');
printf('Parameter          Assignment      Adjusted\n');
printf('---------------------------------------------------------------\n');
printf('umax (1/h)         %9.4f      %9.4f\n', 0.57, umax);
printf('Ks   (g/L)         %9.4f      %9.4f   (fixed)\n', 21, Ks);
printf('kd   (1/h)         %9.4f      %9.4f\n', 0.30, kd);
printf('Yxs  (g/g)         %9.4f      %9.4f   (fixed)\n', Yxs, Yxs);
printf('Yps  (g/g)         %9.4f      %9.4f   (fixed)\n', Yps, Yps);
printf('---------------------------------------------------------------\n');
printf('RMSE X (g/L)       %9.2f      %9.2f\n', rX_nom, rX_aj);
printf('RMSE S (g/L)       %9.2f      %9.2f\n', rS_nom, rS_aj);
printf('RMSE P (g/L)       %9.2f      %9.2f\n', rP_nom, rP_aj);
printf('Global R2          %9.4f      %9.4f\n', R2_nom, R2_aj);
printf('===============================================================\n\n');

printf('Point-by-point comparison (adjusted model):\n');
printf('  t(h)    Xsim    Xexp     Ssim    Sexp     Psim    Pexp\n');
printf('%6.1f %7.2f %7.2f %8.2f %7.2f %8.2f %7.2f\n', ...
       [texp y_aj(:,1) Xexp y_aj(:,2) Sexp y_aj(:,3) Pexp]');
printf('\n');

%===============================================================================
%                   SENSITIVITY ANALYSIS FOR PARAMETER Ks
%===============================================================================
% Justifies the decision to keep Ks = 21 g/L: for each fixed Ks, umax is
% reestimated and the resulting R2 is checked.

printf('Sensitivity analysis to Ks (umax reestimated for each Ks):\n');
printf('   Ks (g/L)     umax (1/h)      R2\n');
Ks_test = [5 10 21 30 50 100];
for i = 1:length(Ks_test)
    Ks = Ks_test(i);
    kd = 0;
    p = abs(fminsearch(@(p) funcao_objetivo([p 0], texp, Xexp, Sexp, Pexp, y0), 0.07, options));
    umax = p;
    [~, y_s] = ode15s(@fun_tessier, texp, y0);
    [~, ~, ~, R2_s] = metricas(y_s, Xexp, Sexp, Pexp);
    printf('%10.1f %14.5f %9.5f\n', Ks, umax, R2_s);
end
printf('\n');

% Restore adjusted parameters for the plots
Ks = 21; umax = p_est(1); kd = p_est(2);

%===============================================================================
%                                PLOTS
%===============================================================================
tit = {'Biomass - X', 'Substrate (glucose) - S', 'Lactic acid - P'};
lab = {'X (g/L)', 'S (g/L)', 'P (g/L)'};
dad = [Xexp Sexp Pexp];

for k = 1:3
    figure(k)
    plot(t_aj, y_aj(:,k), 'b-', 'LineWidth', 2.5)
    hold on
    plot(t_nom, y_nom(:,k), 'r--', 'LineWidth', 2)
    plot(texp, dad(:,k), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 9)
    xlabel('\bf \fontname{Times New Roman} Time (h)', 'Fontsize', 24)
    ylabel(['\bf \fontname{Times New Roman} ' lab{k}], 'Fontsize', 24)
    title(['\bf \fontname{Times New Roman} ' tit{k}], 'Fontsize', 26)
    legend({'Adjusted', 'Assignment', 'Experimental'}, 'Fontsize', 18)
    grid on
    set(gca, 'Fontsize', 18, 'fontWeight', 'bold', 'fontname', 'Times New Roman')
end

printf('Estimation completed.\n');
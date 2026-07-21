function [dy] = fun_tessier(t, y)
    % Função que calcula as derivadas do sistema de fermentação
    % Modelo: Tessier com morte celular
    % y(1) = X (biomassa, g/L)
    % y(2) = S (substrato - glicose, g/L)
    % y(3) = P (produto - ácido lático, g/L)

    global Yxs Yps umax Ks kd

    % Extração das variáveis de estado
    X = y(1);  % Concentração de biomassa (g/L)
    S = y(2);  % Concentração de substrato (g/L)
    P = y(3);  % Concentração de produto (g/L)

    % Proteções numéricas (aplicadas ANTES de calcular mu)
    if X < 1e-6
        X = 1e-6;  % Evita biomassa nula/negativa
    end
    if S < 0
        S = 0;     % Substrato não pode ser negativo
    end

    %=====================================================================
    % MODELO CINÉTICO: TESSIER
    %=====================================================================
    % μ = μmax * (1 - exp(-S/Ks))
    % Características:
    % - Cresce monotonicamente com S
    % - Melhor ajuste para bactérias láticas
    % - Menos sensível a altos valores de S (sem inibição)

    mu = umax * (1 - exp(-S / Ks));

    %=====================================================================
    % SISTEMA DE EQUAÇÕES DIFERENCIAIS
    %=====================================================================
    % Balanço de Massa para BIOMASSA
    % dX/dt = (μ - kd) * X
    % crescimento - morte

    dXdt = (mu - kd) * X;

    %=====================================================================
    % Balanço de Massa para SUBSTRATO
    % dS/dt = -(1/Yxs) * μ * X
    % Consumo de substrato para crescimento celular

    dSdt = -(1 / Yxs) * mu * X;

    %=====================================================================
    % Balanço de Massa para PRODUTO
    % dP/dt = (Yps / Yxs) * μ * X
    % Produção de ácido lático associada ao crescimento

    dPdt = (Yps / Yxs) * mu * X;

    %=====================================================================
    % Montagem do vetor de derivadas
    %=====================================================================
    dy = [dXdt; dSdt; dPdt];

end

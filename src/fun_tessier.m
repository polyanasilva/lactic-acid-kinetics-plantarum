function [dy] = fun_tessier(t, y)
    % Function that calculates the derivatives of the fermentation system
    % Model: Tessier with cell death
    % y(1) = X (biomass, g/L)
    % y(2) = S (substrate - glucose, g/L)
    % y(3) = P (product - lactic acid, g/L)

    global Yxs Yps umax Ks kd

    % State variable extraction
    X = y(1);  % Biomass concentration (g/L)
    S = y(2);  % Substrate concentration (g/L)
    P = y(3);  % Product concentration (g/L)

    % Numerical safeguards (applied BEFORE calculating mu)
    if X < 1e-6
        X = 1e-6;  % Prevents zero or negative biomass
    end
    if S < 0
        S = 0;     % Substrate concentration cannot be negative
    end

    %=====================================================================
    % KINETIC MODEL: TESSIER
    %=====================================================================
    % μ = μmax * (1 - exp(-S/Ks))
    % Characteristics:
    % - Increases monotonically with S
    % - Superior fit for lactic acid bacteria
    % - Less sensitive to high S values (no substrate inhibition)

    mu = umax * (1 - exp(-S / Ks));

    %=====================================================================
    % SYSTEM OF DIFFERENTIAL EQUATIONS
    %=====================================================================
    % Mass Balance for BIOMASS
    % dX/dt = (μ - kd) * X
    % growth - cell death

    dXdt = (mu - kd) * X;

    %=====================================================================
    % Mass Balance for SUBSTRATE
    % dS/dt = -(1/Yxs) * μ * X
    % Substrate consumption for cellular growth

    dSdt = -(1 / Yxs) * mu * X;

    %=====================================================================
    % Mass Balance for PRODUCT
    % dP/dt = (Yps / Yxs) * μ * X
    % Lactic acid production associated with growth

    dPdt = (Yps / Yxs) * mu * X;

    %=====================================================================
    % Assembly of the derivative vector
    %=====================================================================
    dy = [dXdt; dSdt; dPdt];

end
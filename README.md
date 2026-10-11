# Kinetic Modeling of Lactic Acid Production by *Lactiplantibacillus plantarum*

Modeling and simulation of a batch bioreactor culture of *Lactiplantibacillus
plantarum* for lactic acid production, using the **Tessier kinetic model**
for the specific growth rate. Implemented in **GNU Octave**.

Assignment for the course *Computational Modeling and Simulation of Biotechnological
Processes* — Faculty of Biotechnology, UFPA (Prof. Bruno Viegas).

📄 **Full Report:** The complete written report and discussion of results are available in the [artigo.pdf](artigo.pdf) file.

## Model

Specific growth rate (Tessier):

$$\mu = \mu_{max}\left(1 - e^{-S/K_S}\right)$$

System of ODEs from mass balances (batch):

$$\frac{dX}{dt} = (\mu - k_d)\,X \qquad
\frac{dS}{dt} = -\frac{\mu}{Y_{X/S}}\,X \qquad
\frac{dP}{dt} = \frac{\mu\,Y_{P/S}}{Y_{X/S}}\,X$$

where `X` = biomass, `S` = substrate (glucose), and `P` = lactic acid (g/L).

## Results

The nominal parameters from the assignment statement deplete the substrate in less than 1 h (R² = −0.33).
Reestimating only the **kinetic** parameters via least squares (Nelder–Mead simplex),
the fit improves to **R² = 0.985**.

| Parameter | Nominal | Adjusted |
|---|---|---|
| μ_max (h⁻¹) | 0.57 | **0.0705** |
| k_d (h⁻¹) | 0.30 | **≈ 0** |
| K_S (g/L) | 21 | 21 (kept constant) |
| Y_X/S (g/g) | 0.04 | 0.04 (kept constant) |
| Y_P/S (g/g) | 0.43 | 0.43 (kept constant) |

**Dynamics of X, S, and P** — adjusted model (solid line) vs. nominal (dashed) vs. experimental (points):

![Concentration dynamics](figures/fig_fit.png)

**Specific growth rate μ(S)** — nominal overestimates growth across the entire measured range:

![Specific growth rate](figures/fig_mu.png)

**Predicted vs. Experimental Parity** (±15% band):

![Parity plot](figures/fig_parity.png)

## Structure

```
lactic-acid-kinetics-plantarum/
├── src/                     
│   ├── fun_tessier.m            # function containing the system of ODEs (Tessier model)
|   ├── main.m                   # simulation with nominal parameters from the assignment
|   ├── parameter_estimation.m   # parameter reestimation + K_S sensitivity analysis
|   └── generate_figures.m       # generates figures for the figures/ folder
├── figures/                     # result plots (.png and .pdf)
├── artigo.pdf                   # full written report of the project
```

## How to Run

Requires [GNU Octave](https://octave.org/) (tested on version 11.3.0).

```bash
cd src
octave main.m                    # simulation with nominal parameters
octave parameter_estimation.m    # parameter fitting and metrics
octave generate_figures.m        # regenerate figures
```

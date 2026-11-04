# ERCP Timing Before Deterioration in Acute Cholangitis Using Bilirubin WBC Creatinine Platelets and Albumin through a Stochastic Optimal Stopping Model

**Author:** Prihantini  
**Research type:** Literature-parameterized mathematical proof-of-concept  
**Area:** Mathematical modeling, computational medicine, stochastic decision analysis, and optimal stopping

## Overview

This repository provides reproducible Python code, numerical outputs, parameter provenance, and figure scripts for a mathematical proof-of-concept analyzing candidate endoscopic retrograde cholangiopancreatography (ERCP) timing in acute cholangitis. The model evaluates the theoretical trade off between continued stabilization and immediate biliary drainage as simulated laboratory values and latent physiological states change over a 72-hour horizon.

**Important:** No patient level data were used to estimate or validate the model. The numerical results are generated from illustrative initial conditions and provisional coefficients. The reported stopping hours are **not clinical recommendations or validated treatment thresholds**.

## Mathematical and computational components

- Eight state nonlinear ordinary differential equation (ODE) model: total bilirubin, white blood cell count, creatinine, platelets, albumin, inflammatory burden, organ dysfunction, and stability reserve.
- Five continuous severity transformations parameters informed by Tokyo Guidelines 2018 reference criteria.
- Composite dimensionless deterioration index based on prespecified weights.
- Four state continuous time Markov chain (CTMC) describing simulated regime occupancy, the CTMC is descriptive and does not directly enter the Bellman decision.
- Competing waiting related and procedure related **theoretical cost intensities**, not empirically estimated patient event hazards.
- Finite horizon backward Bellman recursion, evaluated on hourly snapshots through 72 hours.
- Eight illustrative starting profile scenarios and a bilirubin creatinine stopping time surface.
- Global uncertainty analysis using 600 Latin hypercube parameter samples and partial rank correlation coefficients (PRCC).
- Local parameter perturbation analysis (±1%) and 1,000 stochastic initial condition perturbations.
- Separate reported post drainage linear relaxation eigenvalues.

## Repository layout

```text
.
├── README.md
├── CITATION.cff
├── requirements.txt
├── cholangitis_optimal_stopping.py
├── make_figures_cholangitis_optimal_stopping_latest.py
├── manuscript/
│   └── Acute Cholangitis Optimal Stopping Article Fix.pdf
└── results/                         # generated after running the model
    ├── Reference_Trajectory.csv
    ├── Reference_Optimal_Stopping_Path.csv
    ├── Scenario_Results.csv
    ├── Bilirubin_Creatinine_Stopping_Surface.csv
    ├── LHS_Parameter_Samples.csv
    ├── Global_Sensitivity_PRCC.csv
    ├── Stochastic_Perturbation_Results.csv
    ├── Local_Elasticity.csv
    ├── Model_Parameters.csv
    ├── Literature_Parameter_Sources.csv
    ├── Post_ERCP_Equilibrium.csv
    ├── Post_ERCP_Eigenvalues.csv
    ├── run_output.json
    ├── figures/                      # created by main script
    └── figures_separate/             # created by figure-only script
```

The `results/` directory is created automatically. If you already have the original published CSV files, retain a separate copy before rerunning the model so you can compare outputs. The manuscript file name in the tree is illustrative, use the actual filename uploaded to the repository.

## Installation

Use Python 3.10 or newer. From the repository root:

```bash
python -m venv .venv
```

Activate the environment:

```bash
# macOS / Linux
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1
```

Install dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Reproduce the numerical analysis

```bash
python cholangitis_optimal_stopping.py
```

The script generates the numerical CSV outputs, `results/run_output.json`, and figure files under `results/figures/`. It also creates DOCX research materials as part of the original project workflow. Depending on your computer, the 600 LHS simulations, 1,000 perturbation simulations, and the bilirubin creatinine sweep may take some time.

The reference ODE uses SciPy `solve_ivp` with RK45, relative tolerance `1e-8`, absolute tolerance `1e-10`, and a 0.25-hour output grid. The repeated uncertainty simulations use a separate fixed-step RK4 routine with a one-hour step and explicit numerical clipping. The Bellman recursion uses hourly decision points. The random seed is `20260819`.

**Reproducibility qualification:** The main script performs an explicit search over `stop_fixed` to select a reference stopping time within its specified target interval. This is numerical tuning against a timing window, **not calibration to patient outcomes**. The seed, numerical methods, software versions, and tuned parameter should be reported when comparing reruns.

## Recreate figures without rerunning simulations

Once the numerical CSV files exist in `results/`, run:

```bash
python make_figures_cholangitis_optimal_stopping_latest.py
```

This creates 14 figures under `results/figures_separate/`. It reads the previously generated CSV files and **does not rerun** the ODE, Latin hypercube sampling, PRCC calculations, or stochastic simulations. Therefore, changing the figure script changes presentation, not the underlying numerical output. The original project ZIP also contains a separately revised figure script, use that script if you specifically want to reproduce its revised figure styling.

## Check key outputs

After the main run, inspect:

```bash
python -c "import json; from pathlib import Path; d=json.loads(Path('results/run_output.json').read_text()); print('Reference stopping hour:',d['reference_stopping_hour']); print('Stochastic:',d['stochastic']); print('Tuned stop_fixed:',d['stop_fixed_tuned'])"
```

Reference results supplied with the project:

| Output | Reported model result |
| --- | ---: |
| Reference candidate stopping time | 29 hours |
| Deterioration index at 24 hours | 0.471 |
| Stability reserve at 24 hours | 0.220 |
| Number of illustrative scenarios | 8 |
| Bilirubin–creatinine input combinations | 483 |
| Latin hypercube parameter samples | 600 |
| Stochastic initial-condition simulations | 1,000 |
| Stochastic median stopping time | 28 hours |
| Stochastic 2.5th–97.5th percentile interval | 19–40 hours |

Small numerical differences may occur with different Python, NumPy, SciPy, or platform versions. Compare both numerical outputs and configuration, rather than only the figures. The 29-hour value is a model specific result, not an empirically supported intervention target.

## Evidence and parameter provenance

`results/Literature_Parameter_Sources.csv` records the literature used to establish clinical reference criteria and contextual ERCP timing evidence. `results/Model_Parameters.csv` distinguishes literature parameters from illustrative initial values and provisional structural or decision coefficients. In particular, published associations and clinical thresholds were **not** converted directly into patient fitted ODE rates or event hazard estimates.

Key parameter context includes Tokyo Guidelines 2018 (Kiriyama et al., 2018, Miura et al., 2018, Mukai et al., 2018, Mayumi et al., 2018), ASGE guidance (Buxbaum et al., 2021), and published timing reviews (Du et al., 2020; Lyu et al., 2022). Full source details and links are recorded in the provenance CSV.

## Interpretation and limitations

This is a simulation mathematical proof-of-concept, not a patient derived prediction model. There is no patient level training set, no prospective validation, no estimated clinical treatment effect, and no demonstrated clinical utility. The latent states, cost intensities, and stopping thresholds depend on model assumptions. PRCC coefficients describe sensitivity within the sampled computational parameter ranges, not causal clinical effects. The CTMC probabilities are model occupancies rather than observed patient transition frequencies. The eigenvalue analysis requires an explicitly defined post drainage vector field and a verified equilibrium before a rigorous local stability claim can be made. The findings should not be used to delay clinically indicated biliary drainage.

## Citation

Use GitHub's **Cite this repository** control when `CITATION.cff` is present. Cite the manuscript separately if it receives a verified publication or DOI. Do not treat the presence of code or a repository citation as evidence of peer reviewed clinical validation.

## Reuse

No license is implied by this README. Add a `LICENSE` file only after choosing the intended reuse terms.

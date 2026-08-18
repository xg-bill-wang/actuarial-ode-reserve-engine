# ODE Model with Life Contingency Application

A Python portfolio project that connects numerical ordinary differential equation methods with a life actuarial reserve application. The model validates Euler, Heun, and fourth-order Runge-Kutta solvers against closed-form ODE examples, then applies the same solver interface to a simplified continuous life insurance reserve model.

The actuarial application uses Thiele's differential equation to show how mortality, interest, premium, and benefit assumptions affect a policy reserve path over time.

## Project Report

The full technical report is available here:

[ODE Model with Life Contingency Application Report](docs/ODE_Model_with_Life_Contingency_Application_Report.pdf)

## What This Project Demonstrates

- Numerical ODE solver implementation in Python.
- Validation against theoretical solutions before actuarial use.
- Life reserve modeling using a simplified Thiele differential equation.
- Scenario testing for mortality, interest-rate, and benefit assumptions.
- Public-facing reporting with charts, tables, and a reproducible script.

## Actuarial Model

The reserve application is based on the following simplified continuous reserve equation:

```text
dV/dt = delta * V + P - mu * (B - V)
```

Where:

- `V` is the policy reserve.
- `delta` is the force of interest.
- `P` is the annual premium rate.
- `mu` is the mortality intensity.
- `B` is the death benefit.

In actuarial terms, the model links reserve growth to investment accumulation, premium inflow, and expected claim pressure. It is intentionally simplified for demonstration and is not intended to replace a full valuation, pricing, or regulatory reserving model.

## Model Workflow

```text
Closed-form ODE examples
        |
        v
Validate Euler, Heun, and RK4 solvers
        |
        v
Apply RK4 solver to Thiele reserve equation
        |
        v
Run assumption sensitivity scenarios
        |
        v
Export CSV outputs and SVG figures
```

## Scenario Outputs

The model compares a base reserve path with three assumption stresses: higher mortality, lower interest, and higher benefit amount.

![Life insurance reserve scenarios](figures/life_reserve_scenarios.svg)

The result is meant to be read actuarially: higher mortality and higher benefit assumptions increase claim pressure, while lower interest reduces reserve accumulation.

## Solver Validation

Before applying the solver to the reserve model, the code tests each numerical method against exponential decay and logistic growth equations with known analytic solutions.

![Solver validation against closed-form ODE solutions](figures/solver_validation.svg)

This validation step shows that the numerical engine behaves as expected before it is used for the actuarial scenario analysis.

## Repository Structure

```text
.
|-- life_reserve_model.py
|-- figures/
|   |-- life_reserve_scenarios.svg
|   `-- solver_validation.svg
|-- docs/
|   |-- ODE_Model_with_Life_Contingency_Application_Report.pdf
|   |-- ODE_Model_with_Life_Contingency_Application_Report.tex
|   `-- ODE_Model_with_Life_Contingency_Application_Report_source.txt
|-- requirements.txt
`-- README.md
```

## How to Run

```bash
python3 life_reserve_model.py
```

Running the script regenerates the model outputs and chart figures.

The current implementation uses only the Python standard library. No external package installation is required for the model script itself.

## Skills Demonstrated

- Python modeling and clean script organization.
- Numerical methods: Euler, Heun, and RK4.
- Actuarial reserve logic and life contingency interpretation.
- Sensitivity testing across actuarial and financial assumptions.
- Model validation against known theoretical results.
- Technical communication for a public GitHub portfolio.

## Limitations and Next Steps

This project is a focused demonstration rather than a production actuarial system. Current limitations include:

- Constant mortality and interest assumptions.
- Simplified premium and benefit structure.
- No expenses, lapses, taxes, capital requirements, or regulatory valuation basis.
- No age-varying life table yet.

The next extension should add a life contingency module with age-based mortality rates, survival probabilities, present value calculations, and deterministic reserves by issue age.

## Credit

This public portfolio project extends ideas from a collaborative ODE engine originally developed with Wenxuan Ding. This repository focuses on the actuarial application layer and presents an independent, recruiter-friendly implementation for public review.



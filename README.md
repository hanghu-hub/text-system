# Power Flow Solvability Conditions: Numerical Validation

This repository provides the MATLAB codes for reproducing the numerical results presented in the paper.

The test cases are based on the standard MATPOWER benchmark systems. All simulations can be reproduced by running the corresponding MATLAB scripts.

## Repository Structure

```
.
├── Table1_SpectralRadius.m
├── Table2_IEEE33.m
├── Table2_IEEE69.m
├── Table2_IEEE85.m
├── Table2_IEEE118.m
├── Table2_IEEE141.m
└── README.md
```

## Numerical Experiments

### Table I: Power-Flow Solutions and Their Spectral-Radius Properties

The script

```
Table1_SpectralRadius.m
```

reproduces the results presented in **Table I**:

> Power-Flow Solutions and Their Spectral-Radius Properties

It calculates multiple power-flow solutions and evaluates their corresponding spectral-radius properties to verify the uniqueness of the spectrally characterized solution.

---

### Table II: Comparison of Certified Loading Factors Under Different Solvability Conditions

The following MATLAB scripts reproduce the results presented in **Table II**:

> Comparison of Certified Loading Factors Under Different Solvability Conditions

Included test systems:

| File | Test System |
|-----|-------------|
| `Table2_IEEE33.m` | IEEE 33-bus radial distribution system |
| `Table2_IEEE69.m` | IEEE 69-bus radial distribution system |
| `Table2_IEEE85.m` | IEEE 85-bus test system |
| `Table2_IEEE118.m` | IEEE 118-bus transmission system |
| `Table2_IEEE141.m` | IEEE 141-bus test system |

These scripts calculate the certified loading factors under different solvability conditions and compare them with the CPF-based solvability boundary.

---

## Figure Reproduction

### Fig. 1: Spectral-radius variation with loading factor

The result of **Fig. 1**:

> Spectral-radius variation with loading factor

is reproduced using:

```
Table2_IEEE33.m
```

The script calculates the variation of the spectral radius with respect to the loading factor and generates the corresponding curve.

---

## Requirements

- MATLAB
- MATPOWER toolbox

The numerical test cases used in this repository are obtained from standard MATPOWER benchmark systems.

---

## How to Run

1. Install MATLAB and configure MATPOWER.

2. Add this repository folder and the MATPOWER folder to the MATLAB path.

3. Run the corresponding MATLAB scripts:

```matlab
Table1_SpectralRadius
```

or

```matlab
Table2_IEEE33
Table2_IEEE69
Table2_IEEE85
Table2_IEEE118
Table2_IEEE141
```

The simulation results will be generated automatically.

---

## Data Source

The numerical test systems used in this repository are obtained from standard MATPOWER benchmark cases.

MATPOWER:
https://matpower.org/

---

## Citation

If you use these codes in your research, please cite the corresponding paper.

```bibtex
@article{
  title={},
  author={},
  journal={},
  year={}
}
```

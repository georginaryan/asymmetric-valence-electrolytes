# Asymmetric Valence Binary Electrolytes
Code accompanying the paper: **"Intermediate-Current Transitions in Asymmetric-Valence Binary Electrolytes"**  
Code author: Georgina Ryan, contact georgina.ryan@maths.ox.ac.uk  
Affiliation: University of Oxford  
Date: 22/04/2026
<img width="191" height="20" alt="image" src="https://github.com/user-attachments/assets/f4ea4bfd-80d0-4707-aa23-bf5ec993a3fb" />

## Overview
This repository contains the analytical and numerical code used to generate the figures in the accompanying paper.

The project combines:
- **Python** for numerical simulation of the time-dependent Poisson–Nernst–Planck (PNP) system  
- **Mathematica** for asymptotic calculations (involving both analytical and numerical methods) and figure generation  
- **MATLAB** for production of Figure 7

All figures from the paper can be reproduced with this code. 

## Dependencies
**Versions used:** Python 3.13.9, MATLAB R2024b, Wolfram Mathematica 14.3.0.0. 
Figures were combined and annotated using Canva.

**Additional packages:**  
Python: numpy, scipy  

## Code Structure and Usage
### Numerical solving (Python)
`PNP_time_evolution.py`  
Python script that solves the 1D Poisson–Nernst–Planck system using a finite-difference discretisation and implicit time integration.  

`generate_time_data.py`  
Driver file for generating numerical solution data using `PNP_time_evolution.py`.  Outputs concentration and potential profiles as CSV files in the `data/` directory, which are used to make Figures 2 and 3.

### Mathematica scripts 
`AnalyticValence.wl`  
Package defining analytical expressions and asymptotic solutions for the cation and anion concentrations and electric potential for the r=1/2,1,2 cases. 

`AsymmetricValence.wl`  
Package defining the numerical method to solve the asymptotic model in the boundary layers for general asymmetric valence electrolytes. 

`FigureSettings.wl`  
Contains common plotting styles, colour schemes, and formatting used across figures.

### Figure generation (Mathematica and MATLAB)
All figures are exported to the `figures/` directory. Final figures were combined and annotated using Canva.

`Figure2.nb`  
Generates Figure 2 from numerical solution data (CSV files).

`Figure3.nb`  
Uses numerical data (CSVs) and functions from AnalyticValence.wl to generate Figure 3.

`Figure4.nb`  
Generates figures using r=2 function from AnalyticValence.wl as the weighted current imbalance J varies.

`Figure5.nb`  
Generates figures using r=2 function from AnalyticValence.wl as the weighted total current I varies.

`Figure6.nb`  
Generates figures using r=1/2,1,2 functions from AnalyticValence.wl and r=1/3,3 from AsymmetricValence.wl.

`Figure7.m`  
Generates the phase diagram of the weighted total current I versus the weighted current imbalance J.

## Running the Code
To reproduce all figures from scratch, first generate the numerical data by running `python scripts/generate_time_data.py` from the project root, which produces the CSV files in the `data/` directory. 
Next, open the Mathematica notebooks in the `notebooks/` directory (`Figure2.nb`–`Figure6.nb`) and evaluate them. 
These notebooks load .wl packages from `source/mathematica/` and numerical data (CSV files) from `data/`, and export their figures to the `figures/` directory. 
Finally, run the MATLAB script `scripts/Figure7.m` to generate the phase diagram, which is also exported to the `figures/` directory. 

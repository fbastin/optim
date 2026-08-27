# Optimization Examples in Julia

This repository stores implementations and application examples of various optimization techniques, primarily written in the **Julia** programming language using **Jupyter Notebooks**.

## Repository Structure

The project is divided into three main categories of optimization:

### 1. LP (Linear Programming)
Contains demonstrations and implementations of foundational linear programming techniques:
- Basic and advanced Simplex methods (`simplexe_basique.ipynb`, `LU et simplexe.ipynb`)
- Affine scaling method (`affineprimal.ipynb`)
- Gauss-Jordan elimination
- Duality theory examples
- Various demonstrations (`Demo1.ipynb` to `Demo11.ipynb`)

### 2. NLP (Non-Linear Programming)
A comprehensive collection of non-linear optimization algorithms and theoretical concepts:
- **Algorithms**: Conjugate Gradient, Newton's method, Steepest Descent, Trust-Region methods (BTR), and Line Search.
- **Applications & Studies**: Deep dive into the Rosenbrock function with multiple visualizations (GIFs and PNGs) showing the iterates of different algorithms.
- **Concepts**: Karush-Kuhn-Tucker (KKT) conditions, Penalty methods, Projected Gradient, and Simulated Annealing.

### 3. SP (Stochastic Programming)
Focuses on optimization under uncertainty:
- The classic Farmer's problem (`Farmer.ipynb`)
- Portfolio chance-constrained programming
- Stochastic Dual Dynamic Programming (SDDP) for hydro management
- Simulation examples (e.g., call center, portfolio)

## Getting Started

To run these notebooks, you will need to have [Julia](https://julialang.org/) installed along with Jupyter. 

The project root contains `Project.toml` and `Manifest.toml` files for Julia's package manager to ensure you have the correct dependencies.

An introductory notebook to the Julia language is also provided (`introJulia.ipynb`).

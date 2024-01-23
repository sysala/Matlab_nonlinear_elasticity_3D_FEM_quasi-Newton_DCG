# Repository for MATLAB Codes acompanying research on Quasi-Newton Iterative Solutions for Non-linear Elasto-Plastic Problems

- Codes for different papers are available in different branches.
---

# Current paper/version: "Quasi-Newton Iterative Solutions for Nonsmooth Elasto-Plastic Problems"

This repository contains MATLAB implementations for our upcoming paper focused on quasi-Newton iterative solutions for nonsmooth elasto-plastic problems. Building upon our previous work in nonlinear elasticity, this new research introduces advanced methods in 3D elasto-plastic models.

- **Title:** Quasi-Newton Iterative Solution Approaches for Nonsmooth Elasto-Plastic Problems
- **Authors:** J. Karátson, S. Sysala, M. Béreš
- **Focus:** Advanced quasi-Newton methods for nonsmooth elasto-plastic problems in 3D, including convergence analysis, novel preconditioners, and iterative solvers.

**Key Features:**
- Application of quasi-Newton methods to nonsmooth elasto-plastic problems in 3D.
- Detailed convergence analysis for the proposed methods.
- Comparative study with standard Newton methods through 3D numerical examples.

### Important note:
- AGMG (Algebraic Multigrid) by Y. Notay is used for experiments. Not freely distributable, academic license available upon request at [AGMG Software and Documentation](http://agmg.eu).
- Once AGMG is available in MATLAB path, you can use it via global variable `AGMG_present = 1` (see test scripts `geocomposite_problem.m` and `homogeneous_problem.m`).
---


# Previous papers:

### Quasi-Newton Variable Preconditioning for Nonlinear Elasticity Systems in 3D
- **DOI:** [10.1002/nla.2537](https://doi.org/10.1002/nla.2537)
- **Branch for Previous Work:** [Paper_10.1002/nla.2537](https://github.com/sysala/matlab_nonlinear_elasticity_3D_FEM_quasi-Newton_DCG/tree/Paper_10.1002/nla.2537)

---
### Acknowledgements:
These codes arise from the codes on elasto-plasticity available at https://github.com/matlabfem/matlab_fem_elastoplasticity

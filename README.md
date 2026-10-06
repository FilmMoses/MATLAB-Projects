# Numerical Analysis & Algorithmic Methods in MATLAB

A computational suite of classical numerical algorithms implemented in MATLAB, focusing on numerical linear algebra, root-finding, interpolation, numerical integration, and initial value problems (IVPs) with rigorous error and convergence verification.

---

## 📌 Repository Overview

This repository contains vectorized, modular MATLAB implementations of foundational numerical methods. Each routine emphasizes computational efficiency, matrix pre-allocation, condition number tracking, and numerical stability:
* **Nonlinear Equation Solvers:** Iterative root-finding routines comparing convergence orders, step counts, and tolerance bounds.
* **Numerical Linear Algebra:** Direct and stationary iterative solvers for systems of linear equations ($Ax = b$), including LU decomposition with partial pivoting and spectral radius analysis.
* **Interpolation & Polynomial Approximation:** High-degree polynomial and piecewise spline interpolation evaluating Runge's phenomenon.
* **Numerical Quadrature & Differentiation:** High-order numerical integration schemes analyzing step size ($h$) asymptotic error scaling.
* **Ordinary Differential Equations (ODEs):** Single-step and multi-stage numerical solvers for first-order and coupled systems of differential equations.

---

## 📐 Mathematical Formulations & Algorithms

### 1. Root-Finding & Nonlinear Solvers
* **Newton-Raphson Method:** Quadratic convergence ($p = 2$) for sufficiently smooth functions:
  $$x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}$$
* **Bisection & Secant Methods:** Robust bracketing and derivative-free approximations measuring convergence rates against stopping criteria:
  $$\frac{\vert{}x_{k+1} - x_k\vert{}}{\vert{}x_{k+1}\vert{}} < \text{tol}$$

### 2. Numerical Linear Algebra
* **LU Decomposition with Partial Pivoting:** Direct factorization solving $PA = LU$, where $P$ is a permutation matrix, $L$ is unit lower triangular, and $U$ is upper triangular.
* **Iterative Solvers (Jacobi & Gauss-Seidel):** Matrix splitting methods $A = M - N$ evaluating spectral radius $\rho(M^{-1}N) < 1$ for convergence:
  $$x^{(k+1)} = D^{-1} \left( b - (L + U)x^{(k)} \right)$$

### 3. Numerical Quadrature
* **Composite Simpson's 1/3 Rule:** 4th-order global truncation error $O(h^4)$:
  $$\int_{a}^{b} f(x) \, dx \approx \frac{h}{3} \left[ f(a) + 2\sum_{j=1}^{n/2-1} f(x_{2j}) + 4\sum_{j=1}^{n/2} f(x_{2j-1}) + f(b) \right]$$

### 4. Initial Value Problems (ODEs)
* **Classical 4th-Order Runge-Kutta (RK4):** Multi-stage integration for non-stiff systems $y' = f(t, y)$:
  $$y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$


## Introduction
This is a high-performance continuation/GMRES (C/GMRES) method for nonlinear model predictive control, HPCGMRES.
This is intended to embedded optimization using high-performance linear algebra package [BLASFEO](https://github.com/giaf/blasfeo.git) and inspired by [HPIPM](https://github.com/giaf/hpipm.git).


## Installation
1. Install [BLASFEO](https://github.com/giaf/blasfeo.git) in your environment as introduced in https://github.com/giaf/blasfeo.git.
2. Clone HPCGMRES in arbitrary directory.
3. 


## License
MIT


## References
- [T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)](https://doi.org/10.1016/j.automatica.2003.11.005)
- [C. T. Kelly, Iterative methods for linear and nonlinear equations, Frontiers in Apllied Mathematics, SIAM (1995)](https://doi.org/10.1137/1.9781611970944)
- [G. Frison, D. Kouzoupis, T. Sartor, A. Zanelli, M. Diehl, BLASFEO: Basic Linear Algebra Subroutines For Embedded Optimization, ACM Transactions on Mathematical Software, 2018](https://arxiv.org/abs/1704.02457)
- [G. Frison, BLASFEO](https://github.com/giaf/blasfeo)
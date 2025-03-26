# Symmetric-GPE-Solver
This code implements the algorithms for a 2D radially symmetric GPE detailed in W. Bao and Y. Cai's paper "MATHEMATICAL THEORY AND NUMERICAL METHODS FOR
BOSE-EINSTEIN CONDENSATION" available on arXiv at https://arxiv.org/pdf/1212.5341.pdf. Results of this code have been benchmarked to known analytic solutions of the quantum harmonic oscillator, and the solutions discussed in  https://arxiv.org/pdf/cond-mat/9608135.pdf in order to obtain numerical estimates of error. 

The 2DRadiallySymmetricGPE_Ground code finds the ground state using BEFD and imaginary time propagation. Numerical errors are on the order of 0.1 dr, where dr is the user specified spatial discretization.
The 2DRadiallySymmetricGPE_Dynamics code numerically time propagates the ground state after quench (in a harmonic trap) using TSFD. Numerical errors again scale with the spatial discretization, which dominates over errors controlled by size of the time propagation step. Errors are on the order of dr at early times and close to the origin. This is the largest magnitude of error. 
If you have already calculated a ground state using other means, or have a saved ground state, you can use it for the time propagation directly.

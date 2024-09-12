# Symmetric-GPE-Solver
This code implements the algorithms for a 2D radially symmetric GPE detailed in W. Bao and Y. Cai's paper "MATHEMATICAL THEORY AND NUMERICAL METHODS FOR
BOSE-EINSTEIN CONDENSATION" available on arXiv at https://arxiv.org/pdf/1212.5341.pdf. 

The first part of this code finds the ground state using BEFD and imaginary time propagation. 
The second part numerically time propagates the ground state after quench (in a harmonic trap) using TSFD. 
If you have already calculated a ground state using other means, or have a saved ground state, you can use it for the time propagation directly.

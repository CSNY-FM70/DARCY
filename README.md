# DARCY
2D Darcy-Pressure Problem 
	div(a(x)grad(p(x)) = f(x) 
	              p(x) = g(x)	on left and right boundary
	      n*grad(p(x)) = 0		on uper and lower boundary

Libraries used:

	- Eigen

Persisting Issues:

	- Optimize nodes & triplets search
	- Make the complete Galerkin Matrix Sparse

To be implemented:

	- Proper Field Generator Class -> Use random/stochastic approach from "An Introduction to Computational Stochastic PDEs" by Lord, Powell, and Shardlow
	- Inhouse Conjugate Gradient solver (or at least create abstract class of Linear::Solvers s.t. one can choose the desired iterative solver) 

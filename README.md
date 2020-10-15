# DARCY
Darcy-Pressure Problem  - Models a unilateral Groundwater Flow in a rectangular 2D-Domain. Numerically solved using a Finite Element approach. Soil permeability -a(x)- is currently only a simple Eigen::Random matrix.
	
	div(a(x)grad(p(x)) = f(x) 
	              
		      p(x) = g(x)	on left and right boundary
	      
	      n*grad(p(x)) = 0		on uper and lower boundary

Libraries used:

	- Eigen (More information - https://eigen.tuxfamily.org/ )

Persisting Issues:

	- Optimize nodes & triplets search -> first hard code optimization *then* use OpenMP
	- Make the complete Galerkin Matrix Sparse

To be implemented:

	- Proper Field Generator Class -> Use random/stochastic approach from "An Introduction to Computational Stochastic PDEs" by Lord, Powell, and Shardlow
	- Inhouse Conjugate Gradient solver (or at least create abstract class of Linear::Solvers s.t. one can choose the desired iterative solver) 

	- CMake Build infrastructure (currently local Make files using GNU c++ compiler)

	- Testing units -> Use Google testing suit or Catch2

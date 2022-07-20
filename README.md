# DARCY
Darcy-Pressure Problem  - Models a unilateral Groundwater Flow in a rectangular 2D-Domain. Numerically solved using a Finite Element approach for a generated soil permeability  field a(x).


$$\nabla\dot(a(x)\nabla(p(x)) = f(x)$$ 
	              
p(x) = g(x)	on left and right boundary
	      
$$\eta\dot\nabla(p(x)) = 0$$	on uper and lower boundary

Pressure distribution over 2D domain and permeability field color-mapped on top.

![Pressure4D](/libdarcy_apps/4D_MAP.png)



Libraries used:

	- Eigen (More information - https://eigen.tuxfamily.org/ )

Persisting Issues:

	- Optimize nodes & triplets search -> first hard code optimization *then* use OpenMP
	- Make the complete Galerkin Matrix Sparse

To be implemented:

	- Inhouse Conjugate Gradient solver (or at least create abstract class of Linear::Solvers s.t. one can choose the desired iterative solver) 

	- Complete CMake Build infrastructure (currently local Make files using GNU c++ compiler)

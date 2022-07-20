# DARCY
Darcy-Pressure Problem  - Models a unilateral Groundwater Flow in a rectangular 2D-Domain for a soil permeability  field a(x). 
Let $$\Omega$$ denote the integration domain and $$\partial \Omega$$ denote the domain boundary, which is segmented in left and right boundary $$\Gamma_0$$
and in the upper and lower boundary $$\Gamma_1$$ for the Neumann boundary condition. where $$\eta$$ is the normal vecor w.r.t. the respective boundary.


$$\nabla\cdot(a(x)\nabla p(x)) = f(x) \quad for x\in\Omega$$ 
	              
$$p(x) = g(x) \quad for x\in\Gamma_0$$
	      
$$\eta\cdot\nabla p(x) = 0 \quad for x\in\Gamma_1$$

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

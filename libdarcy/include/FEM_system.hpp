// Software Lab CES SS2020

// Linear System from Galerkin Approximation
#pragma once

#include <Eigen/Dense>
#include "darcy_system.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/SparseCore>

namespace FEM {
/**
 * @brief Galerkin Linear System
 * @details This class provides the discrete version of the problem from a Finite Elements approach to be solved as a linear system of equations.
 * @date 2020
 * @tparam T Base type to be used for computations.
 */
template<typename T>
class System{
private:
	Darcy::System<T>* _dsy;
	const int ns;
	std::vector<T> bound_data;
	typename Darcy::System<T>::MT PM;
	Eigen::SparseMatrix<T> System_Matrix;
	typename Darcy::System<T>::VT p;
	typename Darcy::System<T>::VT System_RHS;
public:
	template<typename> friend class Solver;
	/**
	* @brief System constructor.
	* @param field Random field used to model system permeability.
	* @param dsy An object of a specialization of the general Darcy::System.
	* @param _ns Number of "squares" in discretization in 1D (amounts to gridpoints-1);
	* @note This constructor is not used in the software as it currently exists.
	*/
	System(typename Darcy::System<T>::MT&, typename Darcy::System<T>&, const int& _ns);

	/**
	* @brief System constructor.
	* @param dsy An object of a specialization of the general Darcy::System.
	* @param _ns Number of "squares" in discretization in 1D (amounts to gridpoints-1);
	*/
	System(typename Darcy::System<T>&,const int & _ns);

	/**
	 * @brief Domain partioned into 2K^2 triangles
	 * @param xv Vector containing x-coordinates of all triangles in the mesh.
	 * @param yv Vector containing y-coordinates of all triangles in the mesh.
	 * @param elt2vert Matrix containing vertex information for all triangles in the mesh.
	 */
	void uniform_mesh(std::vector<T>&, std::vector<T>&, std::vector<std::vector<int>>&);

	/**
	 * @brief Sets up the discretized system.
	 * @details Sets up system matrix and right hand side of the linear system approximating the PDE. Since 0-Neumann boundary conditions exist on part of the domain, the linear system represents all gridpoints in the domain except for the Dirichlet boundary.
	 */
	void setup();

	/**
	* @brief Sets up the permeability matrix from the field matrix.
	* @param _field A matrix modelling a random field.
	*/
	void field_to_perm(typename Darcy::System<T>::MT&);

	/**
	* @brief Galerkin System Matrix A.
	* @return A The system matrix of the discreitzed problem. For K gridpoints in one dimension, this is in R^(K*(K-2)xK*(K-2)) since the Dirichlet boundary points are not part of the discretized solving process.
	*/
	Eigen::SparseMatrix<T>& get_A();

	/**
	* @brief Darcy Pressure.
	* @return p A vector containing the values of the Darcy pressure for each gridpoint. For K gridpoints one dimension, this is in R^(K^2).
	*/
	typename Darcy::System<T>::VT& get_p();

	/**
	* @brief Discretized RHS and boundary data.
	* @return b The right-hand side of the discreitzed problem. For K gridpoints in one dimension, this is in R^(K*(K-2))
	*/
	typename Darcy::System<T>::VT& get_b();

	/**
	* @brief Permeability Matrix
	* @return Matrix for Soil Permeability.
	*/
	typename Darcy::System<T>::MT& get_PM();
};

}

#include"../src/FEM_system.cpp"

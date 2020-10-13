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
 * @tparam T Base type to be used for computations.
 * @tparam K Number of gridpoints in discretized system
 */
template<typename T,int K>
class System{
private:
	typename Darcy::System<T>::MT PM;
	Eigen::SparseMatrix<T> System_Matrix;
	typename Darcy::System<T>::VT p;
	typename Darcy::System<T>::VT System_RHS;
	Darcy::System<T>* _dsy;
public:
	std::vector<int> boundary;
	std::vector<int> interior;
	std::vector<T> bound_data;

	/**
	* @brief System constructor.
	* @param field Random field used to model leading coefficient in Darcy's PDE
	* @param dsy An object of a specialization of the general Darcy::System.
	*/
	System(typename Darcy::System<T>::MT&, typename Darcy::System<T>&);

	/**
	 * @brief Domain partioned into 2K^2 triangles
	 * @param xv Vector containing x-coordinates of all triangles in the mesh.
	 * @param yv Vector containing y-coordinates of all triangles in the mesh.
	 * @param elt2vert Matrix containing vertex information for all triangles in the mesh.
	 */
	void uniform_mesh(std::vector<T>&, std::vector<T>&, std::vector<std::vector<int>>&);
	/**
	 * @brief Initializes A,b,const_a,const_f,nodes and DwB.
	 */
	void setup();

	/**
	* @brief Sets up the permeability matrix from the field matrix.
	* @param _field A matrix modelling a random field.
	*/
	void field_to_perm(typename Darcy::System<T>::MT&);

	/**
	* @brief Galerkin System Matrix A.
	* @return A A matrix which discretizes the continuous problem.
	*/
	Eigen::SparseMatrix<T>& get_A();

	/**
	* @brief Darcy Pressure.
	* @return p A vector containing the values of the Darcy pressure for each gridpoint.
	*/
	typename Darcy::System<T>::VT& get_p();

	/**
	* @brief Discretized RHS and boundary data.
	* @return b A vector containing information regarding the RHS and boundary data for each grid element.
	*/
	typename Darcy::System<T>::VT& get_b();
};

}

#include"../src/FEM_system.cpp"

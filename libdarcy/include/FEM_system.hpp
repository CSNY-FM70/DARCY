// Software Lab CES SS2020
// Gabriel Heredia, Jegor Kravchenko, Wei Yen Kwan, Adrian Lipow
// Supervisor: Sebastian Krumscheid

// Linear System from Galerkin Approximation
#pragma once

#include <Eigen/Dense>
#include "darcy_system.hpp"
#include <cmath>
#include <iostream>

namespace FEM {
/**
 * @brief Galerkin System
 * @details This class provides the discrete version of the problem from a Finite Elements approach to be solved as a linear system of equations.
 * @tparam T Base type for solution precision
 * @tparam K Number of gridpoints in discretized system
 */
template<typename T,int K>
class System{
private:
	typename Darcy::System<T>::MT PM;
	typename Darcy::System<T>::MT A;
	typename Darcy::System<T>::VT p;
	typename Darcy::System<T>::VT b;
	typename Darcy::System<T>::VT const_a;
	typename Darcy::System<T>::VT const_f;
	Darcy::System<T>* _dsy;
public:
	Eigen::Matrix<int,Eigen::Dynamic,1> Dirichlet_nodes;
	Eigen::Matrix<int,Eigen::Dynamic,1> revised_int_nodes;
	typename Darcy::System<T>::VT DwB;
	/**
	* @brief System Constructor
	* @param field Random field to model Equations leading coefficient a
	* @param dsy An object of a specialized Darcy::System 
	*/
	System(typename Darcy::System<T>::MT&, typename Darcy::System<T>&);

	/**
	 * @brief Domain partioned into 2K^2 triangles
	 * @param xv Vector containing discretized x-dimension
	 * @param yv Vector containing discretized y-dimension
	 * @param elt2vert Matrix containing vertex information for all elements
	 */
	void uniform_mesh(typename Darcy::System<T>::VT&, typename Darcy::System<T>::VT&, Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>&);

	/**
	 * @brief Initializes A,b.const_a,const_f,nodes and DwB
	 */
	void setup();

	/**
	* @brief Sets up the permeability matrix from the field matrix
	* @param _field A matrix modelling a random field
	*/
	void field_to_perm(typename Darcy::System<T>::MT&);

	/**
	* @brief Galerkin System Matrix A
	* @return A A matrix which discretizes the continuous problem
	*/
	typename Darcy::System<T>::MT& get_A();

	/**
	* @brief Darcy Pressure
	* @return p A vector containing the values of the darcy pressure for each gridoint after padding in FEM::Solver
	*/
	typename Darcy::System<T>::VT& get_p();

	/**
	* @brief Discretized RHS and boundary data
	* @return b A vector containing information regarding the RHS and boundary data for each grid element
	*/
	typename Darcy::System<T>::VT& get_b();
};

}

#include"../src/FEM_system.cpp"

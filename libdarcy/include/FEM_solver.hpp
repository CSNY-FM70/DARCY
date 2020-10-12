#include "FEM_system.hpp"
#include "linear_system.hpp"
#include "linear_solver.hpp"
#include <Eigen/Dense>

#pragma once

namespace FEM {
/**
* @brief FEM System Solver
* @details This class aims to solve the FEM System using a linear solver from Linear::Solver class
* @tparam T Base type for solution precision
* @tparam K Number of 1D-gridpoints in discretized system
*/

template<typename T,int K>
class Solver{
private:
	typename Linear::Solver<T,(K-1)*(K+1)>& _lsolver;
public:
	/**
	* @brief Solver Constructor
	* @param lsolver Linear solver 
	*/
	Solver(Linear::Solver<T,(K-1)*(K+1)>&);
	/**
	* @brief solve Solves Galerkin System
	* @param femsys An object of type FEM::System to be solved
	*/
	void solve(System<T,K>&);
};

}

#include"../src/FEM_solver.cpp"

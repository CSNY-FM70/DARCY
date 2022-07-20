// Software Lab CES 2020

#include "FEM_system.hpp"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#pragma once

namespace FEM {

/**
* @brief This class aims to solve the FEM System using a Conjugate Gradient Solver.
* @date 2020
* @tparam T Base type to be used for computations.
*/
template<typename T>
class Solver{
private:
	Eigen::ConjugateGradient<Eigen::SparseMatrix<T>,Eigen::Lower|Eigen::Upper> CG;
public:
	/**
	* @brief Solver Constructor.
	*/
	Solver();

	/**
	* @brief solve Solves discretized system.
	* @param femsys An object of type FEM system to be solved.
	*/
	void solve(System<T>& femsys);
};

}

#include"../src/FEM_solver.cpp"

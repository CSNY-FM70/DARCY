#include <iostream>
#include <cassert>
#include "FEM_system.hpp"
#include "FEM_solver.hpp"
#include "darcy_system_test.hpp"
#include "darcy_system_well.hpp"
#include "writeToCsv.cpp"
#include "writeTo4Ddat.cpp"

int main(int argc, const char* argv[]) {
	// assert(argc==2);
	using T = double;
	const int gridpoints = 100;
	const int nvtx = (gridpoints + 1) * (gridpoints + 1);

	std::cout << "Creating boundaries and Darcy system" << std::endl;
	Darcy_Test::System<T>::VTB boundary;
	boundary << 0, 1, 0, 1;
	//Darcy_Test::System<T> dsys(boundary);
	Darcy_Well::System<T> dsys_well(boundary);
	//std::cout << "Successfully created Darcy system" << std::endl;
	T HI = 2.5;
	T LO = -1.5;
	T range = HI - LO;
	typename Darcy::System<T>::MT perm(gridpoints + 1, gridpoints + 1);
	perm = Darcy::System<T>::MT::Random(gridpoints + 1, gridpoints + 1);
	perm = (perm + Eigen::MatrixXd::Constant(gridpoints + 1, gridpoints + 1, 1.)) * range / 2.;
	perm = (perm + Eigen::MatrixXd::Constant(gridpoints + 1, gridpoints + 1, LO));

		std::cout << "FEM system and FEM solver\n";
	//FEM::System<T,gridpoints> femsys(perm,dsys);
	FEM::System<T, gridpoints> femsys_well(perm, dsys_well);
	FEM::Solver<T, gridpoints> femsol;
	// std::cout << "Successfully FEM system and FEM solver" << std::endl;
		std::cout << "Setting up Galerkin System for "<< 2*gridpoints*gridpoints <<" Elements\n";
	femsys_well.setup();
		std::cout << "Solving Galerkin System with "<< nvtx << " Nodes\n";
	femsol.solve(femsys_well);
	
	//Data for 2D Heatmaps
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp_solution = femsys_well.get_p();
	temp_solution.resize(gridpoints + 1, gridpoints + 1);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solution = temp_solution.transpose();
		std::cout << "Writing Data for 2D Heatmaps\n";
	printMatrixToCsv<T>(solution, "solution_well.csv");
	printMatrixToCsv<T>(femsys_well.get_PM(), "PM_well.csv");
	
	//Data for 3D Heatmaps
	femsys_well.get_PM().resize(nvtx, 1);
	Eigen::Matrix<T, Eigen::Dynamic, 1> Sol = femsys_well.get_p();
	Eigen::Matrix<T, Eigen::Dynamic, 1> PM = femsys_well.get_PM();
		std::cout << "Writing Data for 4D Heatmap\n";
	printToDat<T, gridpoints>(Sol, PM, "4DMapping.dat");
	std::cout << "Success!!!\n";
	return 0;
}
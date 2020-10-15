#include <iostream>
#include <cassert>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include "darcy_system_test.hpp"
#include "FEM_system.hpp"
#include "FEM_solver.hpp"
#include "writeToCsv.cpp"


int main(int argc,char* argv[]) {
	assert(argc == 2 && "Insert Samples");
	const int N = std::stoi(argv[1]);
	using T =double;
	const int K = 10;
	Darcy_Test::System<T>::VTB boundary;
 	boundary << 0,1,0,1;
	typename Darcy::System<T>::MT solution;
	T HI = 4;
	T LO = -4;
	T range = HI - LO;
	for(int i=0;i<N;i++){
		//Field Matrix dummy - Added Range
		typename Darcy::System<T>::MT perm(K+1,K+1); 
		perm = Darcy::System<T>::MT::Random(K+1,K+1);
		perm = (perm + Eigen::MatrixXd::Constant(K+1,K+1,1.))*range/2.;
		perm = (perm + Eigen::MatrixXd::Constant(K+1,K+1,LO));
		std::cout << i << ". Field Matrix \n"<< perm << "\n\n";
		Darcy_Test::System<T> dsys_ns(boundary);
		FEM::System<T,K> fsys(perm, dsys_ns);
		std::cout << i << ". Perm Matrix \n" << fsys.get_PM() << "\n\n";
		FEM::Solver<T,K> fsol;
		//fsys.field_to_perm(perm);
		fsys.setup();
		fsol.solve(fsys);
		typename Darcy::System<T>::MT temp_solution = fsys.get_p();
		temp_solution.resize(K+1,K+1);
		solution = temp_solution.transpose();
		std::cout << i << ". Solution \n" << solution << "\n\n";
	}
	//printMatrixToCsv<T>(solution,"solution.csv");
	//std::cout << fsys.get_A() << "\n\n" << fsys.get_b() << "\n\n";
	std::cout << "Success!!" << std::endl;
	return 0;
}

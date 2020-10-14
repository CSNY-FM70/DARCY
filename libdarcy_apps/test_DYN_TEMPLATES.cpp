#include <iostream>
#include <cassert>
#include <Eigen/Dense>

#include "darcy_system_test.hpp"
#include "FEM_system.hpp"
#include "FEM_solver.hpp"

template<typename T,int x, int to>
struct static_for{
	void operator()(){
		typename Darcy_Test::System<T>::VTB boundary;
		boundary << 0,1,0,1;
		Darcy_Test::System<T> dsys_ns(boundary);
		Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> perm = Eigen::Matrix<T,x+1,x+1>::Constant(1);
		FEM::System<T,x> fsys(perm,dsys_ns);
		FEM::Solver<T,x> fsol;
		fsys.setup();
		fsol.solve(fsys);
		//std::cout << "A_int \n" << fsys.get_A() << "\n\n";
		//std::cout << "Drhs \n" << fsys.get_b() << "\n\n";
		std::cout << "Solution for K= " << x << "\n" << fsys.get_p() << "\n\n";
		
		static_for<T,x+1,to>()();
	}
};

template<typename T, int to>
struct static_for<T,to,to>{
	void operator()(){}
};

int main(int argc,char* argv[]) {
	//assert(argc == 2 && "Insert K");
	using T = double;
	static_for<T,5,8>()();
	return 0;
}

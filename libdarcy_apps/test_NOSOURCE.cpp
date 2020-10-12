#include <iostream>
#include <cassert>
#include <Eigen/Dense>

#include "darcy_system_no_source.hpp"
#include "linear_solver_qr.hpp"
#include "FEM_system.hpp"
#include "FEM_solver.hpp"

int main(int argc,char* argv[]) {
	//assert(argc == 2 && "Insert K");
	using T =double;
	const int K = 4;
	const int gp = K+1;
	Darcy_No_Source::System<T>::VTB boundary;
  	boundary << 0,1,0,1;
	const int ne = 2*gp*gp;
	//Field Matrix dummy
	typename Darcy::System<T>::MT perm(gp,gp); perm = Darcy::System<T>::MT::Ones(gp,gp);
	std::cout <<"Perm \n"<< perm << "\n\n";
	Darcy_No_Source::System<T> dsys_ns(boundary);
	//FEM::System<T,K> fsys(perm, dsys_ns);
	FEM::System<T,K> fsys(dsys_ns);
	Linear::Solver_QR<T,(K-1)*(K+1)> lsol;
	FEM::Solver<T,K> fsol(lsol);
	fsys.field_to_perm(perm);
	fsys.setup();
	fsol.solve(fsys);
	
	std::cout << fsys.get_A() << "\n\n" << fsys.get_b() << "\n\n";
	std::cout << fsys.get_p() << "\n\n";
	std::cout << "Success!!" << std::endl;
	return 0;
}

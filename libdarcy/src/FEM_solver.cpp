namespace FEM
{

template <typename T, int K>
Solver<T,K>::Solver(typename Linear::Solver<T,(K-1)*(K+1)>& lsolver):_lsolver(lsolver){}

template<typename T,int K>
void Solver<T,K>::solve(System<T, K>& femsys) {
 	const int n_int = femsys.revised_int_nodes.size();
 	const int n_bound = femsys.Dirichlet_nodes.size();
 	const int nvtx = (K+1)*(K+1);
 	Linear::System<T,(K-1)*(K+1)> lsys(n_int);
 	lsys.A() = femsys.get_A();
 	lsys.b() = femsys.get_b();
 	_lsolver.solve(lsys);
 	femsys.get_p().resize(nvtx,1);
 	//Adding solution for interior nodes
	for(int i=0;i<n_int;i++){
	  	femsys.get_p()(femsys.revised_int_nodes(i)) = lsys.x()(i);
 	}
	//Correcting solution on Dirichlet boundaries
 	for(int j=0;j<n_bound;j++){
	 	 femsys.get_p()(femsys.Dirichlet_nodes(j)) = femsys.DwB(j);
  	}
}

}

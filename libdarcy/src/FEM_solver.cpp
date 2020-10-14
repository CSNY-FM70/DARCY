namespace FEM{

template <typename T, int K>
Solver<T,K>::Solver(){}

template<typename T,int K>
void Solver<T,K>::solve(System<T, K>& femsys) {
 	const int n_int = (K+1)*(K-1);
 	const int n_bound = (K+1)*2;
 	const int nvtx = (K+1)*(K+1);
 	CG.compute(femsys.get_A());
	Eigen::Matrix<T,Eigen::Dynamic,1> u_int = Eigen::Matrix<T,n_int,1>::Zero();
	u_int = CG.solve(femsys.get_b());
 	femsys.get_p().resize(nvtx,1);
 	//Adding solution for interior nodes
	for(int i=0;i<n_int;i++){
	  	femsys.get_p()(femsys.interior.at(i)) = u_int(i);
 	}
	//Correcting solution on Dirichlet boundaries
 	for(int j=0;j<n_bound;j++){
	 	 femsys.get_p()(femsys.boundary.at(j)) = femsys.bound_data.at(j);
  	}
}

}

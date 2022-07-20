namespace FEM{

template<typename T>
Solver<T>::Solver() {}

template<typename T>
void Solver<T>::solve(System<T> & femsys) {
    const int K = femsys.ns;
    const int n_int = (K+1)*(K-1);
 	const int n_bound = (K+1)*2;
 	const int nvtx = (K+1)*(K+1);
 	CG.compute(femsys.System_Matrix);
	Eigen::Matrix<T,Eigen::Dynamic,1> u_int = Darcy::System<T>::VT::Zero(n_int);
	u_int = CG.solve(femsys.System_RHS);
    femsys.p.resize(nvtx,1);
    //Adding solution for interior nodes
    int c = 0;
    for (size_t i = 0; i < K - 1; i++)
        for (size_t j = 0; j < K + 1; j++) {
            femsys.p(c + (K + 1)) = u_int(i + j * (K - 1));
            c++;
        }
	//Correcting solution on Dirichlet boundaries
    int i = 0;
    for (size_t j = 0; j < n_bound; j+=2) {
        femsys.p(i) = femsys.bound_data.at(j);
        femsys.p(i + n_int + (K + 1)) = femsys.bound_data.at(j + 1);
        i++;
    }
}

}

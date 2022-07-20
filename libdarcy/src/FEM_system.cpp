namespace FEM {
template<typename T>
System<T>::System(typename Darcy::System<T>::MT& field, typename Darcy::System<T>& dsy, const int& gridpoints)
	:_dsy(&dsy),ns(gridpoints-1){
		field_to_perm(field);
}
template<typename T>
System<T>::System(typename Darcy::System<T>& dsy, const int& gridpoints)
	:_dsy(&dsy),ns(gridpoints-1){
}
template<typename T>
void System<T>::uniform_mesh(std::vector<T> & xv, std::vector<T> & yv, std::vector<std::vector<int>> & elt2vert) {
	T h=1.0/ns;

	//Discretize axises such that they follow global indexing system
	const int nvtx = (ns+1)*(ns+1);
	xv.resize(nvtx); yv.resize(nvtx);
	for(int i=0;i<ns+1;i++)
		for(int j=0;j<ns+1;j++){
			xv.at(i*(ns+1)+j) = j*h;
			yv.at(i*(ns+1)+j) = i*h;
		}

	const int ne = 2*ns*ns;

	//Global indexing system
	std::vector<std::vector<int>> vv(ns+1,std::vector<int>(ns+1));
	std::vector<std::vector<int>> v1(ns,std::vector<int>(ns)), v2(ns,std::vector<int>(ns)), v3(ns,std::vector<int>(ns)), v4(ns,std::vector<int>(ns));

	typename std::vector<std::vector<int>>::iterator row;
	typename std::vector<int>::iterator col;

	static int c=1;
	for(size_t i=0;i<ns+1;i++)
		for(size_t j=0;j<ns+1;j++) vv.at(j).at(i) = c++;
	c=1;
	for(size_t i=0;i<ns;i++)
		for(size_t j=0;j<ns;j++){
			v1.at(i).at(j) = vv.at(i).at(j);
			v2.at(i).at(j) = vv.at(i+1).at(j);
			v3.at(i).at(j) = vv.at(i).at(j+1);
			v4.at(i).at(j) = vv.at(i+1).at(j+1);
		}


	//List of elements and coresponding vertices based on global indexing system
	elt2vert.resize(ne);
	for(size_t i=0;i<ne;i++)elt2vert.at(i).resize(3);
	int i=0;
	while(i<ne/2){
		int j=0;
		for(size_t l=0;l<ns;l++)
			for(size_t k=0;k<ns;k++){
				//Type A elemets - bottom triangle
				elt2vert.at(i).at(j) = v1.at(k).at(l);
				elt2vert.at(i).at(j+1) = v2.at(k).at(l);
				elt2vert.at(i).at(j+2) = v3.at(k).at(l);
				//Type B elements - upper triangle
				elt2vert.at(i+(ns*ns)).at(j) = v4.at(k).at(l);
				elt2vert.at(i+(ns*ns)).at(j+1) = v3.at(k).at(l);
				elt2vert.at(i+(ns*ns)).at(j+2) = v2.at(k).at(l);
				i++;
			}
	}
}
template<typename T>
void System<T>::setup(){
	const int nvtx=(ns+1)*(ns+1), ne=2*ns*ns;
	T h=1.0/ns;
	std::vector<T> xv, yv;
	std::vector<std::vector<int>> elt2vert;
	uniform_mesh(xv,yv,elt2vert);

	T detJks=h*h;
	//Derivatives of reference element's basis functions
	std::vector<int> dpsi_ds(3), dpsi_dt(3);
	dpsi_ds.at(0)=-1; dpsi_ds.at(1)=1; dpsi_ds.at(2)=0;
	dpsi_dt.at(0)=-1; dpsi_dt.at(1)=0; dpsi_dt.at(2)=1;

	//Piecewise constant data
	std::vector<T> const_a(ne), const_f(ne);
	PM.resize(nvtx,1);
	//Permeability coefficient for each element from GRF.generate
	for(int i=0;i<ne;i++){
		T sum = PM(elt2vert.at(i).at(0)-1)+PM(elt2vert.at(i).at(1)-1)+PM(elt2vert.at(i).at(2)-1);
		const_a.at(i) = (1.0/3.0)*sum;
	}

	//Source term for each element
	for(int i=0;i<ne;i++){
		typename Darcy::System<T>::VTX x;
		T sum = 0;
		for(int j=0;j<3;j++){
			x(0) = xv.at(elt2vert.at(i).at(j)-1);
			x(1) = yv.at(elt2vert.at(i).at(j)-1);
			sum += _dsy->f(x);
		}
		const_f.at(i) = (1.0/3.0) * sum;
	}PM.resize(ns+1,ns+1);

	//Element arrays for piecewise linear finite elements
	std::vector<std::vector<std::vector<T>>> Aks(ne,std::vector<std::vector<T>>(3,std::vector<T>(3)));
	std::vector<std::vector<T>> bks(ne,std::vector<T>(3));
	std::vector<int> delta(ne);
	for(size_t i=0;i<3;i++)
		for(size_t j=0;j<3;j++)
			for(size_t k=0;k<ne/2;k++){
				//Assuming uniform meshes Jacobi matrixes reduce to J_k=diag(ns^-1,ns^-1)
				//delta == grad(psi_p)(x,y) = (J_k^-1)grad(psi_p)(s,t)
				delta.at(k) = (ns*dpsi_ds.at(i))*(ns*dpsi_ds.at(j)) + (ns*dpsi_dt.at(i))*(ns*dpsi_dt.at(j));
				delta.at(k+ne/2) = (-ns*dpsi_ds.at(i))*(-ns*dpsi_ds.at(j)) + (-ns*dpsi_dt.at(i))*(-ns*dpsi_dt.at(j));
				//Adding permeability
				Aks.at(k).at(i).at(j) = detJks*const_a.at(k)*delta.at(k)/2;
				Aks.at(k+ne/2).at(i).at(j) = detJks*const_a.at(k+ne/2)*delta.at(k+ne/2)/2;
			}

	//RHS for each elements vertices
	for(size_t i=0;i<ne;i++){
		bks.at(i).at(0) += const_f.at(i)*(detJks/6);
		bks.at(i).at(1) += const_f.at(i)*(detJks/6);
		bks.at(i).at(2) += const_f.at(i)*(detJks/6);
	}

	//Building Galerkin System from element arrays using global indexing
	std::vector<std::vector<T>> A(nvtx,std::vector<T>(nvtx));
	std::vector<T> b(nvtx);
	for(size_t i=0;i<3;i++){
		for(size_t j=0;j<3;j++)
			for(size_t k=0;k<ne;k++) A.at(elt2vert.at(k).at(i)-1).at(elt2vert.at(k).at(j)-1) += Aks.at(k).at(i).at(j);
		for(size_t l=0;l<ne;l++) b.at(elt2vert.at(l).at(i)-1) += bks.at(l).at(i);
	}

	//Isolating boundary and interior nodes
	std::vector<int> boundary;
	std::vector<int> interior;
	static int current_index = 0;
	boundary.resize((ns+1)*2);
	while(current_index < (ns+1)*2){
		for(size_t i=0;i<nvtx;i++){
			if(xv.at(i)==0 || xv.at(i)==1){
				boundary.at(current_index) = i;
				current_index++;
			}
		}
	}current_index=0;

	std::vector<int> _interior(nvtx);
	for(size_t i=0;i<nvtx;i++){
		_interior.at(i) = i;
		for(size_t j=0;j<(ns+1)*2;j++){
			if(boundary.at(j)==_interior.at(i)) _interior.at(i) = 0;
		}
	}

	for(size_t i=0;i<nvtx;i++){
		if(_interior.at(i)!=0){
			interior.push_back(_interior.at(i));
			current_index++;
		}
	}current_index = 0;

	//Building Linear System to be solved - A_int*u_int = b_int - AD_ib*DwB
	const int n_int = (ns-1)*(ns+1);
	std::vector<T> b_int(n_int);
	for(size_t i=0;i<n_int;i++)b_int.at(i) = b.at(interior.at(i));
	b.resize(0);

	//Creating Sparse Matrix for Interior Nodes
	std::vector<Eigen::Triplet<T>> tripletList;
	tripletList.reserve(n_int);
	for(size_t i=0;i<n_int;i++)
		for(size_t j=0;j<n_int;j++){
			if(A.at(interior.at(i)).at(interior.at(j))!=0){
				tripletList.push_back(Eigen::Triplet<T>(i,j,A.at(interior.at(i)).at(interior.at(j))));
			}
		}
	System_Matrix.resize(n_int, n_int);
	System_Matrix.setFromTriplets(tripletList.begin(), tripletList.end());

	std::vector<std::vector<T>> A_ib(n_int,std::vector<T>((ns+1)*2));
	for(size_t i=0;i<n_int;i++)
		for(size_t j=0;j<(ns+1)*2;j++)A_ib.at(i).at(j) = A.at(interior.at(i)).at(boundary.at(j));

	//Gathering Dirichlet Data
	bound_data.resize((ns+1)*2);
	for(size_t i=0;i<(ns+1)*2;i++){
		typename Darcy::System<T>::VTX x;
		x(0) = xv.at(boundary.at(i));
		x(1) = yv.at(boundary.at(i));
		bound_data.at(i) = _dsy->pD(x);
	}
	xv.resize(0); yv.resize(0);

	//Building right hand side of Linear System
	System_RHS = Darcy::System<T>::VT::Zero(n_int);
	for(size_t i=0;i<n_int;i++){
		for(size_t k=0;k<(ns+1)*2;k++) System_RHS(i) += A_ib.at(i).at(k)*bound_data.at(k);
		System_RHS(i) = b_int.at(i) - System_RHS(i);
	}
	b_int.resize(0);
}
template<typename T>
void System<T>::field_to_perm(typename Darcy::System<T>::MT & _field){
	PM = (exp(_field.array())).matrix();
	//std::cout << "PERMEABILTY \n" << PM << "\n\n";
}
template<typename T>
Eigen::SparseMatrix<T>& System<T>::get_A(){
	return System_Matrix;
}
template<typename T>
 typename Darcy::System<T>::VT& System<T>::get_p() {
	 return p;
 }
 template<typename T>
typename Darcy::System<T>::VT& System<T>::get_b() {
	return System_RHS;
}
template<typename T>
typename Darcy::System<T>::MT& System<T>::get_PM() {
	return PM;
}

}

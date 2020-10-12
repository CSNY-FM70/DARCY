namespace FEM {

template<typename T,int K>
System<T,K>::System(typename Darcy::System<T>::MT& field, typename Darcy::System<T>& dsy)
	:_dsy(&dsy){
		field_to_perm(field);
}

template<typename T,int K>
System<T,K>::System(typename Darcy::System<T>& dsy)
	:_dsy(&dsy){

}

template<typename T,int K>
void System<T,K>::uniform_mesh(typename Darcy::System<T>::VT& xv, typename Darcy::System<T>::VT& yv, Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>& elt2vert){
	const int ns = K;
	T h=1.0/ns;

	//Discretize axises such that they follow global indexing system
	for(int i=0;i<ns+1;i++)
		for(int j=0;j<ns+1;j++)xv((i*(ns+1))+j)=j*h;
	for(int i=0;i<ns+1;i++)
		for(int j=0;j<ns+1;j++)yv((i*(ns+1))+j)=i*h;

	const int ne = 2*ns*ns;

	//Global indexing system
	Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> vv(ns+1,ns+1);
	Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> v1(ns,ns),v2(ns,ns),v3(ns,ns),v4(ns,ns);
	static int c=1;
	for(int i=0;i<(ns+1);i++)
		for(int j=0;j<(ns+1);j++){
			vv(j,i)=c++;
		}c=1;
	for(int i=0;i<ns;i++)
		for(int j=0;j<ns;j++){
			v1(i,j)=vv(i,j);
			v2(i,j)=vv(i+1,j);
			v3(i,j)=vv(i,j+1);
			v4(i,j)=vv(i+1,j+1);
		}
	vv.resize(0,0); //delete vv;

	//List of elements and coresponding vertices based on global indexing system
	int i=0;
	while(i<ne/2){
		int j=0;
		for(int l=0;l<ns;l++)
			for(int k=0;k<ns;k++){
				//Type A elemets - bottom triangle
				elt2vert(i,j) = v1(k,l);
				elt2vert(i,j+1) = v2(k,l);
				elt2vert(i,j+2) = v3(k,l);
				//Type B elements - upper triangle
				elt2vert(i+(ns*ns),j) = v4(k,l);
				elt2vert(i+(ns*ns),j+1) = v3(k,l);
				elt2vert(i+(ns*ns),j+2) = v2(k,l);
				i++;
			}
	}
	v1.resize(0,0); //delete v1;
	v2.resize(0,0); //delete v2;
	v3.resize(0,0); //delete v3;
	v4.resize(0,0); //delete v4;
}

template<typename T,int K>
void System<T,K>::setup(){
	const int ns = K;
	const int nvtx=(ns+1)*(ns+1), ne=2*ns*ns;
	std::cout << "ns = " << ns << std::endl;
	T h=1.0/ns;
	Eigen::Matrix<T,Eigen::Dynamic,1> xv(nvtx); //= Eigen::Matrix<T,nvtx,1>::Zero();
	Eigen::Matrix<T,Eigen::Dynamic,1> yv(nvtx); //= Eigen::Matrix<T,nvtx,1>::Zero();
	Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> elt2vert(ne,3);// = Eigen::Matrix<int,ne,3>::Zero();
	uniform_mesh(xv,yv,elt2vert);

	T detJks=h*h;
	//Derivatives of reference element's basis functions
	Eigen::Matrix<int,1,3> dpsi_ds, dpsi_dt;
	dpsi_ds(0)=-1; dpsi_ds(1)=1; dpsi_ds(2)=0;
	dpsi_dt(0)=-1; dpsi_dt(1)=0; dpsi_dt(2)=1;

	//Piecewise constant data
	const_a.resize(ne,1); //const_a = Eigen::Matrix<T,ne,1>::Zero();
	const_f.resize(ne,1); //const_f = Eigen::Matrix<T,ne,1>::Zero();
	PM.resize(nvtx,1);
	//Permeability coefficient for each element from GRF.generate
	for(int i=0;i<ne;i++){
		int vertices[3];
		vertices[0] = elt2vert(i,0)-1;
		vertices[1] = elt2vert(i,1)-1;
		vertices[2] = elt2vert(i,2)-1;
		T sum = PM(vertices[0])+PM(vertices[1])+PM(vertices[2]);
		const_a(i) = (1.0/3.0)*sum;
	}
	//Source term for each element
	for(int i=0;i<ne;i++){
		typename Darcy::System<T>::VTX x;
		T sum = 0;
		for(int j=0;j<3;j++){
			x(0) = xv(elt2vert(i,j)-1);
			x(1) = yv(elt2vert(i,j)-1);
			sum += _dsy->f(x);
		}
		const_f(i) = (1.0/3.0) * sum;
	}
	
	//Element arrays for piecewise linear finite elements
	Eigen::Matrix<Eigen::Matrix<T,3,3>,Eigen::Dynamic,Eigen::Dynamic> Aks(ne,1);
	Eigen::Matrix<T,3,3> aks;
	for(int i=0;i<ne;i++)Aks(i) = Eigen::Matrix<T,3,3>::Zero();
	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> bks = Eigen::Matrix<T,ne,3>::Zero();
	Eigen::Matrix<int,Eigen::Dynamic,1> delta = Eigen::Matrix<int,ne,1>::Zero();
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			aks = Eigen::Matrix<T,3,3>::Zero();
			for(int k=0;k<ne/2;k++){
				//Assuming uniform meshes Jacobi matrixes reduce to J_k=diag(ns^-1,ns^-1)
				//delta == grad(psi_p)(x,y) = (J_k^-1)grad(psi_p)(s,t)
				delta(k) = (ns*dpsi_ds(i))*(ns*dpsi_ds(j)) + (ns*dpsi_dt(i))*(ns*dpsi_dt(j));
				delta(k+ne/2) = (-ns*dpsi_ds(i))*(-ns*dpsi_ds(j)) + (-ns*dpsi_dt(i))*(-ns*dpsi_dt(j));
				//Adding permeability data
				aks(i,j) = detJks*const_a(k)*delta(k)/2;
				Aks(k) = Aks(k) + aks;
				aks(i,j) = detJks*const_a(k+ne/2)*delta(k+ne/2)/2;
				Aks(k+ne/2) = Aks(k+ne/2) + aks;
			}
		}delta.resize(0);
	//RHS for each elements vertices
	for(int i=0;i<ne;i++){
		bks(i,0) = bks(i,0) + const_f(i)*(detJks/6);
		bks(i,1) = bks(i,1) + const_f(i)*(detJks/6);
		bks(i,2) = bks(i,2) + const_f(i)*(detJks/6);
	}

	//Building Galerkin System from element arrays using global indexing
	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> A = Eigen::Matrix<T,nvtx,nvtx>::Zero();
	Eigen::Matrix<T,Eigen::Dynamic,1> b = Eigen::Matrix<T,nvtx,1>::Zero();
	for(int i=0;i<3;i++){
		Eigen::Matrix<int,Eigen::Dynamic,1> nrow = elt2vert.col(i);
		for(int j=0;j<3;j++){
			Eigen::Matrix<int,Eigen::Dynamic,1> ncol = elt2vert.col(j);
			for(int k=0;k<ne;k++){
			aks = Aks(k);
		        A(nrow(k)-1,ncol(k)-1) += aks(i,j);
			}
		}
		for(int l=0;l<ne;l++){
			b(nrow(l)-1) += bks(l,i);
		}
	}
	//Aks.resize(0,0);
	bks.resize(0,0); elt2vert.resize(0,0);

	//Isolating boundary and interior nodes
	static int current_index = 0;
	Dirichlet_nodes = Eigen::Matrix<int,(ns+1)*2,1>::Zero();
	while(current_index < (ns+1)*2){
		for(int i=0;i<nvtx;i++){
			if(xv(i)==0 || xv(i)==1){
				Dirichlet_nodes(current_index) = i;
				current_index++;
			}
		}
	}current_index=0;

	Eigen::Matrix<int,Eigen::Dynamic,1> _revised_int_nodes(nvtx);
	for(int i=0;i<nvtx;i++){
		_revised_int_nodes(i) = i;
		for(int j=0;j<(ns+1)*2;j++){
			if(Dirichlet_nodes(j)==_revised_int_nodes(i)) _revised_int_nodes(i) = 0;
		}
	}

	for(int i=0;i<nvtx;i++){
		if(_revised_int_nodes(i)!=0){
			revised_int_nodes.conservativeResize(current_index+1);
			revised_int_nodes(current_index)=_revised_int_nodes(i);
			current_index++;
		}
	}current_index = 0; _revised_int_nodes(0);

	//Building Linear System to be solved - A_int*u_int = b_int - AD_ib*DwB
	const int n_int = (ns-1)*(ns+1);
	Eigen::Matrix<T,Eigen::Dynamic,1> b_int = Eigen::Matrix<T,n_int,1>::Zero();
	for(int i=0;i<n_int;i++)b_int(i) = b(revised_int_nodes(i));
	b.resize(0);

	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> A_int = Eigen::Matrix<T,n_int,n_int>::Zero();
	for(int i=0;i<n_int;i++)
		for(int j=0;j<n_int;j++)A_int(i,j) = A(revised_int_nodes(i),revised_int_nodes(j));

	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> AD_ib = Eigen::Matrix<T,n_int,(ns+1)*2>::Zero();
	for(int i=0;i<n_int;i++)
		for(int j=0;j<(ns+1)*2;j++)AD_ib(i,j) = A(revised_int_nodes(i),Dirichlet_nodes(j));
	A.resize(0,0);

	//Gathering Dirichlet Data
	DwB = Eigen::Matrix<T,(ns+1)*2,1>::Zero();
	for(int i=0;i<(ns+1)*2;i++){
		typename Darcy::System<T>::VTX x;
		x(0) = xv(Dirichlet_nodes(i));
		x(1) = yv(Dirichlet_nodes(i));
		DwB(i) = _dsy->pD(x);
	}
	xv.resize(0); yv.resize(0);

	//Building right hand side of Linear System
	Eigen::Matrix<T,Eigen::Dynamic,1> Drhs = Eigen::Matrix<T,n_int,1>::Zero();
	for(int i=0;i<n_int;i++){
		for(int k=0;k<(ns+1)*2;k++) Drhs(i) += AD_ib(i,k)*DwB(k);
		Drhs(i) = b_int(i) - Drhs(i);
	}
	AD_ib.resize(0,0); b_int.resize(0);
	
	get_A() = A_int;
	get_b() = Drhs;
}

template<typename T,int K>
void System<T,K>::field_to_perm(typename Darcy::System<T>::MT& _field) {
	PM = exp(_field.array());
}

template<typename T, int K>
typename Darcy::System<T>::MT& System<T,K>::get_A(){
	return A;
}

template<typename T,int K>
 typename Darcy::System<T>::VT& System<T,K>::get_p(){
	 return p;
 }

template<typename T, int K>
typename Darcy::System<T>::VT& System<T,K>::get_b(){
	return b;
}

}

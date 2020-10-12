#include<iostream>
#include<cassert>
#include<Eigen/Dense>
int main(int argc, char* argv[]){
	assert(argc==2 && "INSERT K");
	const int ns = std::stoi(argv[1]);
	using T =double;
	T h=1.0/ns;
	//std::cout <<"h= "<< h << std::endl;
	const int n2=(ns+1)*(ns+1);	
	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> xv(n2,1);
	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> yv(n2,1);

	for(int i=0;i<ns+1;i++)
		for(int j=0;j<ns+1;j++)xv((i*(ns+1))+j)=h*j;
	for(int i=0;i<ns+1;i++)
		for(int j=0;j<ns+1;j++)yv((i*(ns+1))+j)=h*i;
	//std::cout << xv << "\n\n" << yv <<"\n\n";
	
	const int ne = 2*ns*ns;
	Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> elt2vert(ne,3);
	for(int i=0;i<ne;i++)
		for(int j=0;j<3;j++)elt2vert(i,j)=0;
	
	Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> vv(ns+1,ns+1);
	Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> v1(ns,ns),v2(ns,ns),v3(ns,ns),v4(ns,ns);
	static int c=1;
	for(int i=0;i<(ns+1);i++)
		for(int j=0;j<(ns+1);j++)vv(j,i)=c++;
	for(int i=0;i<ns;i++)
		for(int j=0;j<ns;j++){
			v1(i,j)=vv(i,j);
			v2(i,j)=vv(i+1,j);
			v3(i,j)=vv(i,j+1);
			v4(i,j)=vv(i+1,j+1);
		}c=0;
		//std::cout << v1 << "\n\n";
		//std::cout << v2 << "\n\n";
		//std::cout << v3 << "\n\n";
		//std::cout << v4 << "\n\n";
		//std::cout << vv << "\n\n";
		int i=0;
		while(i<ne/2){
			int j=0;
			for(int l=0;l<ns;l++)
				for(int k=0;k<ns;k++){
					elt2vert(i,j)=v1(k,l);
					elt2vert(i,j+1)=v2(k,l);
					elt2vert(i,j+2)=v3(k,l);
					
					elt2vert(i+(ns*ns),j)=v4(k,l);
					elt2vert(i+(ns*ns),j+1)=v3(k,l);
					elt2vert(i+(ns*ns),j+2)=v2(k,l);
					i++;
				}
		}
		std::cout << elt2vert << "\n\n";
	return 0;
}	

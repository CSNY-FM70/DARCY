#include<iostream>
#include<cassert>

int main(int argc, char *argv[]){
	assert(argc==2 && "INSERT K");
	const int ns = std::stoi(argv[1]);
	using T =  double;
	T h=1.0/ns;
	const int n2=(ns+1)*(ns+1);
	T *xv = new T[n2];
	T *yv = new T[n2];
	for(int i=0;i<ns+1;i++)
		for(int j=0;j<ns+1;j++)xv[i*(ns+1)+j]=h*j;
	for(int i=0;i<ns+1;i++)
		for(int j=0;j<ns+1;j++)yv[i*(ns+1)+j]=i*h;

	//for(int i=0;i<n2;i++)std::cout << xv[i] <<" "<< yv[i] << std::endl;
	delete []xv;
	delete []yv;
	
	const int ne=2*ns*ns;
	int **elt2vert = new int*[ne];
	for(int i=0;i<ne;i++)elt2vert[i]=new int[3];
	
	static int c=1;
	int **vv = new int*[ns+1];
	for(int i=0;i<ns+1;i++){
		vv[i]=new int[ns+1];
		for(int j=0;j<ns+1;j++){
			vv[i][j]=c++;
			//std::cout << vv[i][j] << " ";
		}//std::cout << std::endl;
	}
	
	int **v1 = new int*[ns];
	int **v2 = new int*[ns];
	int **v3 = new int*[ns];
	int **v4 = new int*[ns];
	for(int i=0;i<ns;i++){
		v1[i]=new int[ns];
		v2[i]=new int[ns];
		v3[i]=new int[ns];
		v4[i]=new int[ns];
	}
	for(int i=0;i<ns;i++){
		for(int j=0;j<ns;j++){
			v1[i][j] = vv[i][j];
			v2[i][j] = vv[i+1][j];
			v3[i][j] = vv[i][j+1];
			v4[i][j] = vv[i+1][j+1];

			//std::cout << v3[i][j]<< " ";
		}//std::cout << std::endl;
	}
	for(int i=0;i<ns+1;i++)delete[]vv[i];
	delete []vv;
	
	int i=0;
	while(i<ne/2){
		int j=0;
		for(int l=0;l<ns;l++)
			for(int k=0;k<ns;k++){
				elt2vert[i][j] = v1[l][k];
				elt2vert[i][j+1] = v3[l][k];
				elt2vert[i][j+2] = v2[l][k];

				elt2vert[i+(ns*ns)][j] = v4[l][k];
				elt2vert[i+(ns*ns)][j+1] = v2[l][k];
				elt2vert[i+(ns*ns)][j+2] = v3[l][k];
				i++;
			}
	}
	for (int i=0;i<ns;i++){
		delete []v1[i];
		delete []v2[i];
		delete []v3[i];
		delete []v4[i];
	}
	delete []v1;
	delete []v2;
	delete []v3;
	delete []v4;

	for(int i=0;i<ne;i++)std::cout << elt2vert[i][0] << "\t" << elt2vert[i][1]<<"\t"<<elt2vert[i][2]<<"\n";
	for(int i=0;i<ne;i++)delete []elt2vert[i];
	delete []elt2vert;
	return 0;
}

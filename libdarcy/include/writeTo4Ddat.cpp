#include <iostream>
#include <fstream>
#include <cmath>
#include "Eigen/Dense"
// This function writes solution and permeability into a 4 columns file.
// Inputs consist of 2 vectors for solution, permeability and a string for filename.
// example: printToDat<T,K>(femsys.get_p(),femsys.get_PM(), "data.dat");
template<typename T,int K>
void printToDat(const Eigen::Matrix<T,Eigen::Dynamic,1>& Sol, const Eigen::Matrix<T,Eigen::Dynamic,1>& PM, const std::string fileName) {
	std::ofstream stream;
        stream.open(fileName);
    	const int nvtx = (K+1)*(K+1);
        const int ns = K;
	T h = 1.0/ns;
        Eigen::Matrix<T,nvtx,1> xv, yv;
        for(int i=0;i<ns+1;i++)
    		for(int j=0;j<ns+1;j++){
        		xv(i*(ns+1)+j) = j*h;
                	yv(i*(ns+1)+j) = i*h;
       		}	

    	for (int i = 0; i < ns+1; i++) {
		for(int j=0; j<ns+1;j++){
        		stream << xv(i*(ns+1)+j) << "\t" << yv(i*(ns+1)+j) << "\t" << Sol(i*(ns+1)+j) << "\t" << PM(i*(ns+1)+j) << "\n";
        	}
		stream << "\n";
	}
}

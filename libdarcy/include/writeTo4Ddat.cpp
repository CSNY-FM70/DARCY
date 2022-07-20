#include <iostream>
#include <fstream>
#include <cmath>
#include "Eigen/Dense"
// This function writes solution and permeability into a 4 columns file.
// Inputs consist of 2 vectors for solution, permeability and a string for filename.
// example: printToDat(femsys.get_p(),femsys.get_PM(), "data.dat");
template<typename T>
void print_to_dat(const int& K,const Eigen::Matrix<T,Eigen::Dynamic,1>& Sol, const Eigen::Matrix<T,Eigen::Dynamic,1>& PM, const std::string fileName) {
	std::ofstream stream;
        stream.open(fileName);
    	const int nvtx = K*K;
        const int ns = K-1;
	T h = 1.0/ns;
        Eigen::VectorXd xv(nvtx), yv(nvtx);
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

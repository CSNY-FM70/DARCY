#include <iostream>
#include <fstream>
#include "Eigen/Dense"
// This function writes an Eigen matrix into a .csv (comma separated values) file.
// Inputs consist of a matrix and a filename as a string.
// example: printMatrixToCsv<T>(matrix, "data.csv");
template<typename T>
void printMatrixToCsv(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M, const std::string fileName) {
    std::ofstream stream;
    stream.open(fileName);
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            stream << M(i,j) << ",";
        }
        stream << "\n";
    }
}

// Library of mathematical functions to produce random fields.

#pragma once

#include "Eigen/Dense"
#include "unsupported/Eigen/FFT"
#include <iostream>
#include <math.h>
#include <random>
#include <complex>
#include "darcy_system.hpp"
#include <vector>

namespace GRF_matlib {

    template<typename T>
    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename T>
    using MTC = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename T>
    using VT = Eigen::Matrix<T, 1, Eigen::Dynamic>;

    template<typename T>
    using VTX=Eigen::Matrix<T,2,1>;


    /**
     * @brief Fast Fourier transformation of a 2D matrix.
     * @date 2020
     * @tparam T Base type to be used for computations.
     * @param A Input matrix.
     * @return Fourier transofermed matrix.
     */
    template<typename T>
    MTC<T> fft(const MTC<T>& A) {

        Eigen::FFT<T> fft;
        MTC<T> B(A.rows(), A.cols());
        std::vector<std::complex<T>> *aux_1 = new std::vector<std::complex<T>>(A.cols());
        std::vector<std::complex<T>> *aux_2 = new std::vector<std::complex<T>>(A.cols());

        for (int k = 0; k < A.rows(); k++) {
            for (int i = 0; i < A.cols(); i++) {
                aux_1->at(i) = A(k,i);
            }
            fft.fwd(*aux_2,*aux_1);
            for (int i = 0; i < A.cols(); i++) {
                B(k,i) = aux_2->at(i);
            }
        }
        delete aux_1;
        delete aux_2;

        aux_1 = new std::vector<std::complex<T>>(A.rows());
        aux_2 = new std::vector<std::complex<T>>(A.rows());

        for (int j = 0; j < A.cols(); j++) {
            for (int i = 0; i < A.rows(); i++) {
                aux_1->at(i) = B(i,j);
            }
            fft.fwd(*aux_2,*aux_1);
            for (int i = 0; i < A.rows(); i++) {
                B(i,j) = aux_2->at(i);
            }
        }
        delete aux_1;
        delete aux_2;

        return B;

    }


    /**
     * @brief Inverse fast Fourier transformation of a 2D matrix.
     * @date 2020
     * @tparam T Base type to be used for computations.
     * @param A Input matrix.
     * @return Inverse Fourier transofermed matrix.
     */
    template<typename T>
    MTC<T> ifft(const MT<T>& A) {

        Eigen::FFT<T> fft;
        MTC<T> B(A.rows(), A.cols());
        std::vector<std::complex<T>> *aux_1 = new std::vector<std::complex<T>>(A.cols());
        std::vector<std::complex<T>> *aux_2 = new std::vector<std::complex<T>>(A.cols());

        for (int k = 0; k < A.rows(); k++) {
            for (int i = 0; i < A.cols(); i++) {
                aux_1->at(i) = A(k,i);
            }
            fft.inv(*aux_2,*aux_1);
            for (int i = 0; i < A.cols(); i++) {
                B(k,i) = aux_2->at(i);
            }
        }
        delete aux_1;
        delete aux_2;

        aux_1 = new std::vector<std::complex<T>>(A.rows());
        aux_2 = new std::vector<std::complex<T>>(A.rows());
        for (int j = 0; j < A.cols(); j++) {
            for (int i = 0; i < A.rows(); i++) {
                aux_1->at(i) = B(i,j);
            }
            fft.inv(*aux_2,*aux_1);
            for (int i = 0; i < A.rows(); i++) {
                B(i,j) = aux_2->at(i);
            }
        }
        delete aux_1;
        delete aux_2;

        return B;

    }


    /** @brief fft_shift shifts zero-frequency component to center of spectrum.
     *  For a matrix it swaps the first quadrant with the third and the second with the forth.
     *
     *     I  II    fft-shift    IV III
     *    III IV   ---------->   II  I
     *
     *  In case of uneven number of rows or columns, quadrant I gets the bigger "half".
     * @date 2020
     * @tparam T Base type to be used for computations.
     * @param A Input matrix.
     * @return Fourier transofermed matrix.
     */    
    template<typename T>
    MT<T> fft_shift(const MT<T>& A) {

        MT<T> M;
        M.resize(A.rows(), A.cols());
        int row_shift, col_shift;
        // Check if number of rows is odd or even, set rwo_shift accordingly.
        if (A.rows()%2 == 0) {
            row_shift = A.rows()/2;
        }
        else {
            row_shift = A.rows()/2+1;
        }
        // Check if number of columns is odd or even, set rwo_shift accordingly.
        if (A.cols()%2 == 0) {
            col_shift = A.cols()/2;
        }
        else {
            col_shift = A.cols()/2+1;
        }

        // write IV into I
        for (int i = 0; i < A.rows()-row_shift; i++) {
            for (int j = 0; j < A.cols()-col_shift; j++) {
                M(i,j) = A(i+row_shift, j+col_shift);
            }
        }
        // write I into IV
        for (int i = A.rows()-row_shift; i < A.rows(); i++) {
            for (int j = A.cols()-col_shift; j < A.cols(); j++) {
                M(i,j) = A(i-A.rows()+row_shift, j-A.cols()+col_shift);
            }
        }
        // write III into II
        for (int i = A.rows()-row_shift; i < A.rows(); i++) {
            for (int j = 0; j < A.cols()-col_shift; j++) {
                M(i,j) = A(i-A.rows()+row_shift, j+col_shift);
            }
        }
        // write II into III
        for (int i = 0; i < A.rows()-row_shift; i++) {
            for (int j = A.cols()-col_shift; j < A.cols(); j++) {
                M(i,j) = A(i+row_shift, j-A.cols()+col_shift);
            }
        }

        return M;

    }

    /**
     * @brief Assemble reduced covariance matrix, based on covariance function provided by related darcy system.
     * @date 2020
     * @tparam T Base type to be used for computations.
     * @param n_1 Number of rows of reduced covariance matrix.
     * @param n_2 Number of columns of reduced covariance matrix.
     * @param dx_1 Discretization size in vertical space.
     * @param dx_2 Discretization size in horizontal space.
     * @param d Reference to related darcy system.
     * @return Reduced covariance matrix.
     */
    template<typename T>
    MT<T> reduced_cov(const int n_1, const int n_2, const T dx_1, const T dx_2, Darcy::System<T>* d) {

        MT<T> c_red(2*n_1-1,2*n_2-1);
        VTX<T> eval_point;

        for (int i = 0; i < 2*n_1-1; i++) {
            eval_point(0) = (i-n_1+1)*dx_1;
            for (int j = 0; j < 2*n_2-1; j++) {
                eval_point(1) = (j-n_2+1)*dx_2;
                c_red(i,j) = d->c(eval_point);
            }
        }
        return c_red;

    }

    /**
     * @brief Generate random sample.
     * @date 2020
     * @tparam T Base type to be used for computations.
     * @param c_red Reduced covariance matrix.
     * @param n_1 Number of rows of the sample.
     * @param n_2 Number of columns of the sample.
     * @param seed Seed for random generator.
     * @return A random sample.
     */
    template<typename T>
    MT<T> circ_cov_sample_2d(const MT<T>& c_red, const int n_1, const int n_2, int seed) {

        // random generator with normal distribution
        std::mt19937 generator;    
        generator.seed(seed);        
        std::normal_distribution<T> distribution(0,1);

        // c_red is of size n_1 x n_2
        int N = n_1*n_2;
        MTC<T> lambda = N*ifft(c_red);

        //generate normal distributed random numbers
        for (int i = 0; i < n_1; i++) {
            for (int j = 0; j < n_2; j++) {
                lambda(i,j) = sqrt(lambda(i,j)) * distribution(generator);
            }
        }
        lambda = (1.0/sqrt(N))*fft(lambda);
        lambda.resize(lambda.rows()*lambda.cols(),1);
        MT<T> lambda_real = lambda.array().real();
        lambda.resize(0,0);
        return lambda_real;

    }

    /**
     * @brief Extract field from random samples.
     * @date 2020
     * @tparam T Base type to be used for computations.
     * @param c_red Reduced covariance matrix.
     * @param n_1 Number of rows of the sample.
     * @param n_2 Number of columns of the sample.
     * @param m_1 Padding parameter for vertical dimension.
     * @param m_2 Padding parameter for horizontal dimension.
     * @param seed Seed for random generator.
     * @return A random field.
     */
    template<typename T>
    MT<T> circ_embeded_sample_2dB(const MT<T>& c_red, const int n_1, const int n_2, const int m_1, const int m_2, int seed) {

        const int nn_1 = n_1 + m_1;
        const int nn_2 = n_2 + m_2;
        // form reduced matrix of BCCB extension of BTTB matrix C*
        MT<T> tilde_c_red, aux;
        // tilde_c_red gets a zero row and zero column on the top left side.
        tilde_c_red.resize(2*nn_1, 2*nn_2);

        for (int i = 1; i < tilde_c_red.rows(); i++) {
            for (int j = 1; j < tilde_c_red.cols(); j++) {
                tilde_c_red(i,j) = c_red(i-1,j-1);
            }
        }        
        for (int i = 0; i < tilde_c_red.rows(); i++) {
            tilde_c_red(i,0) = 0;
        }        
        for (int j = 1; j < tilde_c_red.cols(); j++) {
            tilde_c_red(0,j) = 0;
        }        
        
        tilde_c_red=fft_shift(tilde_c_red);
        MT<T> u = circ_cov_sample_2d(tilde_c_red, 2*nn_1, 2*nn_2, seed);
        // put columns of u into a vector
        u.resize(4*nn_1*nn_2,1);
        // discard entries past 2*nn_1*n_2
        aux = u;
        u.resize(2*nn_1*n_2,1);
        for (int i = 0; i < 2*nn_1*n_2; i++) {
            u(i) = aux(i);
        }
        u.resize(nn_1, 2*n_2); // resize u
        // use only the first n_1 rows and every second column
        aux = u;
        if (u.cols()%2 == 0) u.resize(n_1, u.cols()/2);
        else u.resize(n_1, u.cols()/2+1);
        for (int i = 0; i < u.rows(); i++) {
            for (int j = 0; j < u.cols(); j++) {
                u(i,j) = aux(i,j*2);
            }
        }
        return u;

    }

}

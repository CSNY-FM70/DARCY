// Software Lab CES SS2020

#pragma once

#include "Eigen/Dense"
#include "GRF.hpp"
#include "FEM_system.hpp"
#include "FEM_solver.hpp"
#include "darcy_system.hpp"

namespace UQ {
  /**
   * @brief Class used to asses stochastic effects of uncertainty.
   * @details This class performs a Monte Carlo sample of the discretized Darcy system calculates mean pressure at each gridpoint as well as the pressure variance at each gridpoint.
   * @date 2020
   * @tparam T Base type to be used for computations.
   */
  template<typename T>
  class Sampler {
  private:
    //using MT = typename Darcy::System<T>::MT;
    using VT = typename Darcy::System<T>::VT;

    int N;
    typename Darcy::System<T>::MT mean, sum_square_difference;
    GRF::Generator<T>* generator;
    int iterations;
    FEM::System<T>* system;
    FEM::Solver<T>* solver;
    const int gridpoints;

    /**
     * @brief Performs a single Monte Carlo evaluation.
     * @return Returns the solved Darcy System.
     */
    typename Darcy::System<T>::MT single_sample();

    /**
     * @brief Updates both the mean and variance over all iterations so far.
     * @param current_solution The solution matrix to be added to the sample.
     */
    void update_mean_and_variance(const typename Darcy::System<T>::MT& current_solution);

  public:
    /**
     * @brief Sampler Constructor for an object of type Sampler with a given number of Monte Carlo iterations.
     * @param N Number of Monte Carlo samples to be performed.
     * @param sys Reference to discretized system to be solved.
     * @param sol Reference to solver for discretized system.
     * @param _K Number of gridpoints in 1 dimension.
     */
    Sampler(int N, GRF::Generator<T>& gen, FEM::System<T>& sys, FEM::Solver<T>& sol, const int& _K);

    /**
     * @brief sample Runs a Monte Carlo sample.
     */
    void sample();

    /**
     * @brief get_mean The mean pressure of the system at each gridpoint.
     * @return A matrix containing the mean pressure at each gridpoint of the system.
     */
    const typename Darcy::System<T>::MT& get_mean();

    /**
     * @brief get_variance The pressure variance of the system at each gridpoint.
     * @return A matrix containing the pressure variance at each gridpoint of the system.
     */
    const typename Darcy::System<T>::MT get_variance();

    /**
     * @brief get_N Number of Monte Carlo iterations to be performed.
     * @return The currently set number of Monte Carlo iterations.
     */
    int get_N();

    /**
     * @brief set_N Sets the number of Monte Carlo iterations to be performed.
     * @param new_N Number of iterations.
     */
    void set_N(const int new_N);

    /**
    * @brief Reset everything in the sampler to 0.
    */
    void reset();
  };

}

#include "../src/sampler.cpp"

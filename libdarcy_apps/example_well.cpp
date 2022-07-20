//----------------------------------------------------
//CES Software Lab 2020
//Example program 2
//This program models a Darcy system that contains an extracting well.
//Command line inputs are: Number of Monte Carlo samples, sigma (models area affected by well), jet (models strength of well), l (covariance parameter) and number of gridpoints in one dimension
//The random generator here is seeded with a constant, so that this program is deterministic and can be run as a test case when the software is first installed.
//Output files are CSV format: permeability field of the last sampled system, mean pressure of the system, pressure variance of the system.
//----------------------------------------------------

#include <iostream>
#include <cassert>
#include "darcy_system_well.hpp"
#include "Darcy.hpp"


int main(int argc, const char* argv[]) {
  assert(argc==6 && "Insert #samples,sigma(0.01-2),jet(-100 - 100) and l(0.1-0.9)->cov(r,l),Gridpoints(10-200)");
  try {
    using T = double;
    if (atoi(argv[1]) < 1)
      throw std::runtime_error("number of samples must be positive!");
    int samples = atoi(argv[1]);
  // 1. Create Darcy system boundaries and then the system
  // 2. Create the linear solver, FEM system and FEM solver
  // 3. Create the generator
  // 4. Create the sampler
  // 5. Run Monte Carlo sample
  // 6. Return all the info

    Darcy_Well::System<T>::VTB boundary;
    boundary << 0,1,0,1;
    if (atof(argv[2]) < 0.01 || atof(argv[2]) > 2)
      throw std::runtime_error("sigma must be between 0.01 and 2!");
    T sigma = static_cast<T>(atof(argv[2]));
    if (atof(argv[3]) < -100 || atof(argv[3]) > 100)
      throw std::runtime_error("jet must be between -100 and 100!");
    T jet = static_cast<T>(atof(argv[3]));
    if (atof(argv[4]) < 0.1 || atof(argv[4]) > 0.9)
      throw std::runtime_error("covariance parameter l must be between 0.1 and 0.9!");
    T lcov = static_cast<T>(atof(argv[4]));
    if (atoi(argv[5]) < 2 || atoi(argv[5]) > 300)
        throw std::runtime_error("No less than 2 Gridpoint and if considerably more than 200 check available memory");
    const int gridpoints = atoi(argv[5]);
    Darcy_Well::System<T> dsys_well(boundary,sigma,jet,lcov);
    FEM::System<T> femsys_well(dsys_well,gridpoints);
    FEM::Solver<T> femsol;
    bool deterministic = true;
    GRF::Generator<T> generator_well(gridpoints, dsys_well, deterministic);
    UQ::Sampler<T> sampler_well(samples, generator_well, femsys_well, femsol, gridpoints);
    std::cout << "Starting Monte Carlo sample" << std::endl;
    sampler_well.sample();
    std::cout << "Successfully finished Monte Carlo sampling" << std::endl;

    print_matrix_to_csv<T>(femsys_well.get_PM(), "PM_well.csv");
    print_matrix_to_csv<T>(sampler_well.get_mean(), "mean_well.csv");
    print_matrix_to_csv<T>(sampler_well.get_variance(), "variance_well.csv");
  } catch (std::exception& ex) {
    std::cout << ex.what() << std::endl;
  }
  return 0;
}

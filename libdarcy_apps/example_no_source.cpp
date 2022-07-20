//----------------------------------------------------
//CES Software Lab 2020
//Example program 1
//This program models a Darcy system that has no source.
//Command line inputs are: Number of Monte Carlo samples and number of gridpoints in one dimension
//The random generator here is seeded with a constant, so that this program is deterministic and can be run as a test case when the software is first installed.
//Output files are CSV format: permeability field of the last sampled system, mean pressure of the system, pressure variance of the system.
//----------------------------------------------------

#include <iostream>
#include <cassert>
#include "Darcy.hpp"
#include "darcy_system_no_source.hpp"


int main(int argc, const char* argv[]) {
  assert(argc==4 && "Insert #samples, l(0.1-0.9)->cov(r,l), Gridpoints(10-200)");
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

    Darcy_No_Source::System<T>::VTB boundary;
    boundary << 0,1,0,1;
    if (atof(argv[2]) < 0.1 || atof(argv[2]) > 0.9)
      throw std::runtime_error("covariance parameter l must be between 0.1 and 0.9!");
    T lcov = static_cast<T>(atof(argv[2]));
    if (atoi(argv[3]) < 2 || atoi(argv[3]) > 300)
        throw std::runtime_error("No less than 2 Gridpoint and if considerably more than 200 check available memory");
    const int gridpoints = atoi(argv[3]);
    Darcy_No_Source::System<T> dsys(boundary,lcov);
    FEM::System<T> femsys(dsys,gridpoints);
    FEM::Solver<T> femsol;
    bool deterministic = true;
    GRF::Generator<T> generator(gridpoints, dsys, deterministic);
    UQ::Sampler<T> sampler(samples, generator, femsys, femsol, gridpoints);
    std::cout << "Starting Monte Carlo sample" << std::endl;
    sampler.sample();
    std::cout << "Successfully finished Monte Carlo sampling" << std::endl;

    print_matrix_to_csv(femsys.get_PM(), "PM_no_source.csv");
    print_matrix_to_csv(sampler.get_mean(), "mean_no_source.csv");
    print_matrix_to_csv(sampler.get_variance(), "variance_no_source.csv");
  } catch (std::exception& ex) {
    std::cout << ex.what() << std::endl;
  }
  return 0;
}

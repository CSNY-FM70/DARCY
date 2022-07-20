
#include <iostream>
#include <cassert>
#include "darcy_system_well.hpp"
#include "Darcy.hpp"

int main(int argc, const char* argv[]) {
	// assert(argc==2);
	assert(argc == 6 && "Insert #samples,sigma(0.01-2),jet(-100 - 100) and l(0.1-0.9)->cov(r,l),Gridpoints(10-200)");
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
		boundary << 0, 1, 0, 1;
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
		Darcy_Well::System<T> dsys_well(boundary, sigma, jet, lcov);
		std::cout << "Successfully created Darcy system\n" << std::endl;
		FEM::System<T> femsys_well(dsys_well, gridpoints);
		FEM::Solver<T> femsol;
		std::cout << "FEM system and FEM solver instantieted\n";
		bool deterministic = true;
		GRF::Generator<T> generator_well(gridpoints, dsys_well, deterministic);
		UQ::Sampler<T> sampler_well(samples, generator_well, femsys_well, femsol, gridpoints);
		std::cout << "Starting Monte Carlo sample\n" << std::endl;
		sampler_well.sample();
		std::cout << "Successfully finished Monte Carlo sampling\n" << std::endl;

		const int nvtx = (gridpoints + 1.0) * (gridpoints + 1.0);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sol = sampler_well.get_mean();
		Sol.resize(nvtx, 1);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PM = femsys_well.get_PM();
		PM.resize(nvtx, 1);
		//Data for 3D Heatmaps
		std::cout << "Writing Data for 4D Heatmap\n";
		print_to_dat<T>(gridpoints+1,Sol, PM, "4DMapping.dat");
		std::cout << "Success!!!\n";
	}
	catch (std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}
	return 0;
}
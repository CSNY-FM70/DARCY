
#include <iostream>
#include <cassert>
#include "darcy_system_well.hpp"
#include "Darcy.hpp"

int main(int argc, const char* argv[]) {
	// assert(argc==2);
	assert(argc == 5 && "Insert sigma(0.01-2),jet(-100 - 100) and l(0.1-0.9)->cov(r,l),Gridpoints(10-200)");
	try {
		using T = double;
		// 1. Create Darcy system boundaries and then the system
		// 2. Create the linear solver, FEM system and FEM solver
		// 3. Create the generator
		// 4. Create the sampler
		// 5. Run Monte Carlo sample
		// 6. Return all the info
		Darcy_Well::System<T>::VTB boundary;
		boundary << 0, 1, 0, 1;
		if (atof(argv[1]) < 0.01 || atof(argv[1]) > 2)
			throw std::runtime_error("sigma must be between 0.01 and 2!");
		T sigma = static_cast<T>(atof(argv[1]));
		if (atof(argv[2]) < -100 || atof(argv[2]) > 100)
			throw std::runtime_error("jet must be between -100 and 100!");
		T jet = static_cast<T>(atof(argv[2]));
		if (atof(argv[3]) < 0.1 || atof(argv[3]) > 0.9)
			throw std::runtime_error("covariance parameter l must be between 0.1 and 0.9!");
		T lcov = static_cast<T>(atof(argv[3]));
		if (atoi(argv[4]) < 2 || atoi(argv[4]) > 300)
			throw std::runtime_error("No less than 2 Gridpoint and if considerably more than 200 check available memory");
		const int gridpoints = atoi(argv[4]);
		Darcy_Well::System<T> dsys_well(boundary, sigma, jet, lcov);
		std::cout << "Successfully created Darcy system\n" << std::endl;
		FEM::System<T> femsys_well(dsys_well, gridpoints);
		FEM::Solver<T> femsol;
		std::cout << "FEM system and FEM solver instantieted\n";
		bool deterministic = true;
		GRF::Generator<T> generator_well(gridpoints, dsys_well, deterministic);
		
		generator_well.generate_field();
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> random_sample = generator_well.get_field();
		femsys_well.field_to_perm(random_sample);
		femsys_well.setup();
		femsol.solve(femsys_well);

		const int nvtx = gridpoints*gridpoints;

		Eigen::Matrix<T, Eigen::Dynamic, 1> Sol = femsys_well.get_p();
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> PM = femsys_well.get_PM();
		PM.resize(nvtx, 1);

		//Data for 3D Heatmaps
		std::cout << "Writing Data for 4D Heatmap\n";
		print_to_dat<T>(gridpoints,Sol, PM, "4DMapping.dat");
		std::cout << "Success!!!\n";
	}
	catch (std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}
	return 0;
}
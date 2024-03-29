// Software Lab CES SS2020

// header file for random field generator

#pragma once

#include "Eigen/Dense"
#include "darcy_system.hpp"
#include <chrono>

namespace GRF {

    /**
    * @brief Generator This class generates Gaussian random fields.
    * @date 2020
    * @tparam T Base type to be used for computations.
    */
    template<typename T>
    class Generator {

        private:

        using MT=Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
        using VTX=Eigen::Matrix<T,2,1>;

        int gridpoints;
        MT field;
        Darcy::System<T>* dsys;
        bool deterministic;
        int seed;

        public :

        /**
         * @brief Generator Constructor for Generator objects.
         * @param p_gridpoints The number of gridpoints in 1 dimension.
         * @param d A reference to a Darcy system, required for getting the covariance function of the system.
         * @param det can be set to true for getting reproducible results, is set to false by default.
         */
        Generator(const int p_gridpoints, Darcy::System<T>& d, bool det = false);

        /**
         * @brief generate_field Generates a Gaussian random field.
         */
        void generate_field();

        /**
         * @brief get_field Returns the last field generated by the generator.
         * @return The last field generated.
         */
        MT& get_field();

    };

}

#include "../src/GRF.cpp"

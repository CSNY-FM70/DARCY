// Example case with with a Gaussain function as source and constant pressure on left and right boundaries

#include "darcy_system.hpp"
#include <cmath>
//#include "Eigen/Dense"

namespace Darcy_Well {
/**
   * @brief Test case modeling an extracting well.
   * @details Sample implementation of a Darcy system. Uses constant boundary pressure and a RHS that models an extracting well.
   * @tparam T Base type to be used for computations.
   */
  template<typename T>
  class System : public Darcy::System<T> {
    using Darcy::System<T>::boundary;
  public:

    /**
     * @brief System Constructor
     * @param bound Vector containing the boundary locations (rectangular domain)
     */
    System(typename Darcy::System<T>::VTB bound);

    /**
     * @brief f Right-hand side of the system. Models an extracting well centered at the middle of the domain.
     * @param x Point at which f is to be evaluated.
     * @return Value of RHS at given point.
     */
    T f(typename Darcy::System<T>::VTX x);

    /**
     * @brief pD Dirichlet boundary condition of the system. Is constant along left and right boundary.
     * @param x Point at which boundary condition is to be evaluated.
     * @return Pressure at given point. -1 if x is not on the boundary.
     */
    T pD(typename Darcy::System<T>::VTX x);
  };
}

#include "../src/darcy_system_well.cpp"

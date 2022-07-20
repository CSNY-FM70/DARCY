// Software Lab CES SS2020

// Example case with no sources and constant pressure on left and right boundaries

#include "darcy_system.hpp"
#include <cmath>
//#include "Eigen/Dense"

namespace Darcy_Well {
/**
   * @brief Test case modeling an extracting well.
   * @details Sample implementation of a Darcy system. Uses constant boundary pressure and a RHS that models an extracting well.
   * @date 2020
   * @tparam T Base type to be used for computations.
   */
  template<typename T>
  class System : public Darcy::System<T> {
    
    using Darcy::System<T>::boundary;
    const T sigma;
    const T jet;
    const T l;

  public:

    /**
     * @brief System Constructor
     * @param bound Vector containing the boundary locations (rectangular domain)
     * @param _sigma parameter to model area affected by well.
     * @param _jet parameter to model strength of the well.
     * @param _l covariance parameter.
     */
    System(const typename Darcy::System<T>::VTB&, const T& _sigma,const T& _jet,const T& _l);

    /**
     * @brief f Right-hand side of the system. Models an extracting well centered at the middle of the domain.
     * @param x Point at which f is to be evaluated.
     * @return Value of RHS at given point.
     */
    T f(const typename Darcy::System<T>::VTX& x);

    /**
     * @brief pD Dirichlet boundary condition of the system. Is constant along left and right boundary.
     * @param x Point at which boundary condition is to be evaluated.
     * @return Pressure at given point. -1 if x is not on the boundary.
     */
    T pD(const typename Darcy::System<T>::VTX& x);

    /**
     * @brief c Covariance function of the system.
     * @param dist Directional vector between two points.
     * @return Covariance for given points.
     */
    T c(const typename Darcy::System<T>::VTX& dist);
  };
}

#include "../src/darcy_system_well.cpp"

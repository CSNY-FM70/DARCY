// Software Lab CES SS2020
// Abstract base class from which the user-specified RHS and boundary conditions
// of the Darcy system are to be specified.
#pragma once

#include <Eigen/Dense>

namespace Darcy {
  /**
   * @brief Abstract base class from which a Darcy system can be specialized.
   * @details This class provides the basic interface through which a Darcy system is implemented. The system is defined by the right-hand side, the Dirichlet boundary conditions and the isotropic covariance.
   * @date 2020
   * @tparam T Base type to be used for computations.
   */
  template<typename T>
  class System {
  public:

    using VT=Eigen::Matrix<T,Eigen::Dynamic,1>;
    using VTX=Eigen::Matrix<T,2,1>;
    using VTB=Eigen::Matrix<T,4,1>;
    using MT=Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>; //TODO: is this necessary?

  protected:
    VT boundary;

  public:
    /**
     * @brief System Constructor
     * @param bound Vector containing the boundary locations (rectangular domain), format should be xL,xU,yL,yU
     */
    System(const VTB& bound);

    /**
     * @brief f Right-hand side of the system.
     * @param x Point at which f is to be evaluated.
     * @return Value of RHS at given point.
     * @note Remeber to add "typename Darcy::System<T>" to VTX in parameter when specializing the class.
     */
    virtual T f(const VTX& x) = 0;

    /**
     * @brief pD Dirichlet boundary condition of the system.
     * @param x Point at which boundary condition is to be evaluated.
     * @return Pressure at given point.
     * @note Implementation should handle case where x is not on the boundary.
     * @note Remeber to add "typename Darcy::System<T>" to VTX in parameter when specializing the class.
     */
    virtual T pD(const VTX& x) = 0;

    /**
     * @brief c Covariance function of the system.
     * @param dist Directional vector between two points.
     * @return Covariance for given points.
     * @note Remeber to add "typename Darcy::System<T>" to VTX in parameter when specializing the class.
     */
    virtual T c(const VTX& dist) = 0;

  };
}

#include "../src/darcy_system.cpp"

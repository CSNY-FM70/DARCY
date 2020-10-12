#pragma once

#include <Eigen/Dense>

namespace Darcy {
  /**
   * @brief Abstract base class from which a Darcy system can be specialized.
   * @details This class provides the basic interface through which a Darcy system is implemented. The system is defined by the right-hand side, the Dirichlet boundary conditions(Neumann boundaries homogeneous)
   * @tparam T Base type for solution precision
   */
template<typename T>
class System {
public:
    using VT=Eigen::Matrix<T,Eigen::Dynamic,1>;
    using VTX=Eigen::Matrix<T,2,1>;
    using VTB=Eigen::Matrix<T,4,1>;
    using MT=Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

protected:
    VT boundary;

public:
    /**
     * @brief System Constructor
     * @param bound Vector containing the boundary locations (rectangular domain), format should be xL,xU,yL,yU
     */
    System(VTB bound);

    /**
     * @brief f Right-hand side of the system.
     * @param x Point at which f is to be evaluated.
     * @return Value of RHS at given point.
     */
    virtual T f(VTX x) = 0;

    /**
     * @brief pD Dirichlet boundary condition of the system.
     * @param x Point at which boundary condition is to be evaluated.
     * @return Pressure at given point.
     * @note Implementation should handle case where x is not on the Dirichlet boundaries.
     */
    virtual T pD(VTX x) = 0;
  };
}

#include "../src/darcy_system.cpp"

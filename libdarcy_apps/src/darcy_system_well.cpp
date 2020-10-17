namespace Darcy_Well {

template<typename T>
System<T>::System(typename Darcy::System<T>::VTB bound) : Darcy::System<T>(bound){}

template<typename T>
T System<T>::f(typename Darcy::System<T>::VTX x) {
  T x_middle = 0.5;
  T y_middle = 0.5;
  T x_part = pow(x(0)-x_middle,2)/0.05;
  T y_part = pow(x(1)-y_middle,2)/0.05;
  return -100*exp(-(x_part+y_part));
}

template<typename T>
T System<T>::pD(typename Darcy::System<T>::VTX x) {
  if (x(0) == boundary(0)) return 5; // Lower x bound
  else if (x(0) == boundary(1)) return 4; // Upper x bound
  // Negative value to indicate error if given point is not on the left or right boundary
  else return -1;
}

}

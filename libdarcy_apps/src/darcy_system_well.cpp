namespace Darcy_Well {

template<typename T>
System<T>::System(const typename Darcy::System<T>::VTB& bound, const T& _sigma, const T& _jet, const T& _l) : Darcy::System<T>(bound), sigma(_sigma), jet(_jet),l(_l){}

template<typename T>
T System<T>::f(const typename Darcy::System<T>::VTX& x) {
  T x_middle = 0.5;
  T y_middle = 0.5;
  T x_part = pow(x(0)-x_middle,2)/sigma;
  T y_part = pow(x(1)-y_middle,2)/sigma;
  return jet*exp(-(x_part+y_part));
}

template<typename T>
T System<T>::pD(const typename Darcy::System<T>::VTX& x) {
  if (x(0) == boundary(0)) return 0; // Lower x bound
  else if (x(0) == boundary(1)) return 0; // Upper x bound
  else return -1;
}

template<typename T>
T System<T>::c(const typename Darcy::System<T>::VTX& dist) {
  T r = dist.norm();
  T r_square = pow(r,2);
  T l_square = pow(l,2);
  //std::cout << "c(" << dist << ")=" << exp(-rsquare/lsquare) << std::endl;
  return exp(-r_square/l_square);
  //T result = exp(-rsquare/lsquare);
  //if (result < pow(10,-10)) result = 0;
  //return result;
}

}

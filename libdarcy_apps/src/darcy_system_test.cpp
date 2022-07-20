namespace Darcy_Test {

template<typename T>
System<T>::System(const typename Darcy::System<T>::VTB& bound) : Darcy::System<T>(bound){}

template<typename T>
T System<T>::f(const typename Darcy::System<T>::VTX& x) {
  return 0;
}

template<typename T>
T System<T>::pD(const typename Darcy::System<T>::VTX& x) {
  if (x(0) == boundary(0)) return 10; // Lower x bound
  else if (x(0) == boundary(1)) return 0; // Upper x bound
  else return -1;
}

template<typename T>
T System<T>::c(const typename Darcy::System<T>::VTX& dist) {
  T r = dist.norm();
  T r_square = pow(r,2);
  T l_square = pow(l,2);
  return exp(-r_square/l_square);
  

}

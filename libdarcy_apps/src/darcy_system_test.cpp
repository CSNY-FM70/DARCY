namespace Darcy_Test {

template<typename T>
System<T>::System(typename Darcy::System<T>::VTB bound) : Darcy::System<T>(bound){}

template<typename T>
T System<T>::f(typename Darcy::System<T>::VTX x) {
  return 1;
}

template<typename T>
T System<T>::pD(typename Darcy::System<T>::VTX x) {
  if (x(0) == boundary(0)) return 0; // Lower x bound
  else if (x(0) == boundary(1)) return 0; // Upper x bound
  else return -1;
}

}

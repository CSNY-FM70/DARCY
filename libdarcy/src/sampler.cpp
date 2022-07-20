namespace UQ {

  template<typename T>
  Sampler<T>::Sampler(int N, GRF::Generator<T>& gen, FEM::System<T>& sys, FEM::Solver<T>& sol,const int& _K) : N(N), generator(&gen), iterations(0), system(&sys), solver(&sol),gridpoints(_K){
    mean = Darcy::System<T>::MT::Zero(gridpoints,gridpoints);
    sum_square_difference = Darcy::System<T>::MT::Zero(gridpoints,gridpoints);
  }

  template<typename T>
  typename Darcy::System<T>::MT Sampler<T>::single_sample() {
    generator->generate_field();
    typename Darcy::System<T>::MT field = generator->get_field();
    system->field_to_perm(field);
    system->setup();
    solver->solve(*system);
    typename Darcy::System<T>::MT solution = system->get_p();
    solution.resize(gridpoints, gridpoints);
    return solution;
  }

  template<typename T>
  void Sampler<T>::update_mean_and_variance(const typename Darcy::System<T>::MT& current_solution) {
    iterations++;
    typename Darcy::System<T>::MT delta = current_solution - mean;
    mean += delta/iterations;
    typename Darcy::System<T>::MT delta2 = current_solution - mean;
    sum_square_difference += (delta.array()*delta2.array()).matrix();
  }

  template<typename T>
  void Sampler<T>::sample() {
    for (int i = 0; i < N; i++) {
      typename Darcy::System<T>::MT solved_system = single_sample();
      update_mean_and_variance(solved_system);
    }
  }

  template<typename T>
  const typename Darcy::System<T>::MT& Sampler<T>::get_mean() {
    return mean;
  }

  template<typename T>
  const typename Darcy::System<T>::MT Sampler<T>::get_variance() {
    return sum_square_difference/iterations;
  }

  template<typename T>
  int Sampler<T>::get_N() {
    return N;
  }

  template<typename T>
  void Sampler<T>::set_N(const int new_N) {
    N = new_N;
  }

  template<typename T>
  void Sampler<T>::reset() {
    N = 0;
    iterations = 0;
    mean = typename Darcy::System<T>::MT::Zero(gridpoints,gridpoints);
    sum_square_difference = typename Darcy::System<T>::MT::Zero(gridpoints,gridpoints);
  }
}

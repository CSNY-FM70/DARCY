#include "GRF_matlib.cpp"

namespace GRF {

    template<typename T>
    Generator<T>::Generator(const int p_gridpoints, Darcy::System<T>& d, bool det) : gridpoints(p_gridpoints), dsys(&d), deterministic(det), seed(-1) {}

    template<typename T>
    void Generator<T>::generate_field() {
        if (deterministic) seed++;
        else seed = std::chrono::system_clock::now().time_since_epoch().count();
        const int m = 8*gridpoints; //padding parameter
        const T dx = 1.0/(gridpoints-1);
        MT reduced_covariance = GRF_matlib::reduced_cov<T>(gridpoints+m, gridpoints+m, dx, dx, dsys);
        field = GRF_matlib::circ_embeded_sample_2dB<T>(reduced_covariance, gridpoints, gridpoints, m, m, seed);
    }

    template<typename T>
    typename Generator<T>::MT& Generator<T>::get_field() {
        return field;
    }

}

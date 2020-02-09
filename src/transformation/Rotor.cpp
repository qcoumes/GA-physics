#include <transformation/Rotor.hpp>


namespace transformation {
    
    template<typename T>
    Rotor<T>::Rotor(const c3ga::Mvec<T> &bivector, T angle) :
        Versor<T>(
            cos(0.5 * angle) - bivector * sin(0.5 * angle),
            c3ga::Mvec<T>(cos(0.5 * angle) - bivector * sin(0.5 * angle)).inv()
        ) {
    }
}

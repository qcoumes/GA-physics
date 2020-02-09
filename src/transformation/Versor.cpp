#include <transformation/Versor.hpp>


namespace transformation {
    
    template<typename T>
    Versor<T>::Versor(const c3ga::Mvec <T> &t_versor, const c3ga::Mvec <T> &t_inverse) :
        versor(t_versor), inverse(this->translator.inv()) {
    }
    
    
    template<typename T>
    c3ga::Mvec <T> Versor<T>::operator()(const c3ga::Mvec <T> &object) {
        return versor * object * inverse;
    }
}

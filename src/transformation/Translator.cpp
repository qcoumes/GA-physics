#include <transformation/Translator.hpp>


namespace transformation {
    
    template<typename T>
    Translator<T>::Translator(const c3ga::Mvec<T> &vector) :
        Versor<T>(1 - 0.5 * vector * c3ga::ei<T>(), c3ga::Mvec<T>(1 - 0.5 * vector * c3ga::ei<T>()).inv()) {
    }
}

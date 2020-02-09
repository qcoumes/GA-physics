#include <transformation/Dilator.hpp>


namespace transformation {
    
    template<typename T>
    Dilator<T>::Dilator(T factor) :
        Versor<T>(
            1 - (1 - factor) / (1 + factor) * c3ga::e0i<T>(),
            c3ga::Mvec<T>(1 - (1 - factor) / (1 + factor) * c3ga::e0i<T>()).inv()
        ) {
    }
}

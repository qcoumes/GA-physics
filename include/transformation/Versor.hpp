#ifndef GA_PHYSICS_VERSOR_HPP
#define GA_PHYSICS_VERSOR_HPP

#include <c3ga/Mvec.hpp>


namespace transformation {
    
    template<typename T>
    class Versor {
    
        protected:
            c3ga::Mvec<T> versor;
            c3ga::Mvec<T> inverse;
        
        public:
            
            explicit Versor(const c3ga::Mvec<T> &versor, const c3ga::Mvec<T> &inverse);
            
            c3ga::Mvec<T> operator()(const c3ga::Mvec<T> &object);
    };
}

#endif //GA_PHYSICS_VERSOR_HPP

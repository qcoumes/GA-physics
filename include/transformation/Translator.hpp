#ifndef GA_PHYSICS_TRANSLATOR_HPP
#define GA_PHYSICS_TRANSLATOR_HPP

#include <c3ga/Mvec.hpp>

#include <transformation/Versor.hpp>


namespace transformation {
    
    template<typename T>
    class Translator : public Versor<T> {
        
        public:
            
            explicit Translator(const c3ga::Mvec<T> &vector);
    };
    
    
    
    
}

#endif //GA_PHYSICS_TRANSLATOR_HPP

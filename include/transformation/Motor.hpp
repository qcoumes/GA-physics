#ifndef GA_PHYSICS_TRANSLATOR_HPP
#define GA_PHYSICS_TRANSLATOR_HPP

#include <c3ga/Mvec.hpp>

#include <transformation/Versor.hpp>


namespace transformation {
    
    template<typename T>
    class Motor : public Versor<T> {
        
        public:
            
            Motor(const std::initializer_list<T> &versors);
    };
    
    
    
    
}

#endif //GA_PHYSICS_TRANSLATOR_HPP

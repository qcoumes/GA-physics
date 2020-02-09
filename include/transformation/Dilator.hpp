#ifndef GA_PHYSICS_DILATOR_HPP
#define GA_PHYSICS_DILATOR_HPP

#include <c3ga/Mvec.hpp>

#include <transformation/Versor.hpp>


namespace transformation {
    
    template<typename T>
    class Dilator : public Versor<T> {
        
        public:
            
            explicit Dilator(T factor);
    };
    
    
    
    
}

#endif //GA_PHYSICS_DILATOR_HPP

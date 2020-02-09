#ifndef GA_PHYSICS_ROTOR_HPP
#define GA_PHYSICS_ROTOR_HPP

#include <c3ga/Mvec.hpp>

#include <transformation/Versor.hpp>


namespace transformation {
    
    template<typename T>
    class Rotor : public Versor<T> {
        
        public:
            
            explicit Rotor(const c3ga::Mvec<T> &bivector, T angle);
    };
}

#endif // GA_PHYSICS_ROTOR_HPP

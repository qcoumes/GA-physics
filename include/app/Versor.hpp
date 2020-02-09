#ifndef GA_PHYSICS_VERSOR_HPP
#define GA_PHYSICS_VERSOR_HPP

#include <c3ga/Mvec.hpp>


namespace app {
    
    template<typename T>
    class Versor {
        
        private:
            c3ga::Mvec<T> versor;
            c3ga::Mvec<T> inverse;
            
            explicit Versor(const c3ga::Mvec<T> &versor);
            
            Versor(const c3ga::Mvec<T> &versor, const c3ga::Mvec<T> &inverse);
        
        public:
            
            static Versor<T> dilator(T factor);
            
            static Versor<T> rotor(const c3ga::Mvec<T> &bivector, T angle);
            
            static Versor<T> translator(const c3ga::Mvec<T> &vector);
            
            [[nodiscard]] Versor<T> operator*(const c3ga::Mvec<T> &other) const;
            
            [[nodiscard]] c3ga::Mvec<T> operator()(const c3ga::Mvec<T> &object) const;
    };
    
    
    
    template<typename T>
    Versor<T>::Versor(const c3ga::Mvec<T> &t_versor) :
            versor(t_versor), inverse(this->versor.inv()) {
    }
    
    
    template<typename T>
    Versor<T>::Versor(const c3ga::Mvec<T> &t_versor, const c3ga::Mvec<T> &t_inverse) :
            versor(t_versor), inverse(t_inverse) {
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::dilator(T factor) {
        return Versor<T>(1 - (1 - factor) / (1 + factor) * c3ga::e0i<T>());
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::rotor(const c3ga::Mvec<T> &bivector, T angle) {
        return Versor<T>(cos(0.5 * angle) - bivector * sin(0.5 * angle));
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::translator(const c3ga::Mvec<T> &vector) {
        return Versor<T>(1 - 0.5 * vector * c3ga::ei<T>());
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::operator*(const c3ga::Mvec<T> &other) const {
        return Versor<T>(this->versor * other.versor, this->inverse * other.inverse);
    }
    
    
    template<typename T>
    c3ga::Mvec<T> Versor<T>::operator()(const c3ga::Mvec<T> &object) const {
        return versor * object * inverse;
    }
    
    /////////////////////////// DEDUCTION GUIDES ///////////////////////////////
    
    template<class T>
    Versor(const c3ga::Mvec<T> &) -> Versor<T>;
    
    template<class T>
    Versor(const c3ga::Mvec<T> &, const c3ga::Mvec<T> &) -> Versor<T>;
}

#endif //GA_PHYSICS_VERSOR_HPP

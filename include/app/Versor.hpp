#ifndef GA_PHYSICS_VERSOR_HPP
#define GA_PHYSICS_VERSOR_HPP

#include <c3ga/Mvec.hpp>


namespace app {
    
    template<typename T>
    class Versor {
        
        private:
            c3ga::Mvec<T> versor;
            c3ga::Mvec<T> inverse;
        
            explicit Versor(const c3ga::Mvec<T> &mvec);
        
            Versor(const c3ga::Mvec<T> &mvec, const c3ga::Mvec<T> &inverse);
        
        public:
        
            Versor() = default;
            
            static Versor<T> reflector(const c3ga::Mvec<T> &mvec);
            
            static Versor<T> dilator(T factor);
            
            static Versor<T> rotor(const c3ga::Mvec<T> &bivector, double angle);
            
            static Versor<T> translator(const c3ga::Mvec<T> &vector);
            
            static Versor<T> translator(T x, T y, T z);
            
            [[nodiscard]] Versor<T> operator*(const Versor<T> &other) const;
            
            [[nodiscard]] c3ga::Mvec<T> operator()(const c3ga::Mvec<T> &object) const;
            
            [[nodiscard]] Versor<T> operator()(const Versor<T> &other) const;
    };
    
    
    
    template<typename T>
    Versor<T>::Versor(const c3ga::Mvec<T> &mvec) :
            versor(mvec), inverse(this->versor.inv()) {
    }
    
    
    template<typename T>
    Versor<T>::Versor(const c3ga::Mvec<T> &mvec, const c3ga::Mvec<T> &t_inverse) :
            versor(mvec), inverse(t_inverse) {
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::reflector(const c3ga::Mvec<T> &mvec) {
        return Versor<T>(mvec);
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::dilator(T factor) {
        return Versor<T>(1 - (1 - factor) / (1 + factor) * c3ga::e0i<T>());
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::rotor(const c3ga::Mvec<T> &bivector, double angle) {
        return Versor<T>(cos(0.5 * angle) - bivector * sin(0.5 * angle));
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::translator(const c3ga::Mvec<T> &vector) {
        return Versor<T>(1 - 0.5 * vector * c3ga::ei<T>());
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::translator(T x, T y, T z) {
        return Versor<T>(1 - 0.5 * (x * c3ga::e1<T>() + y * c3ga::e2<T>() + z * c3ga::e3<T>()) * c3ga::ei<T>());
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::operator*(const Versor<T> &other) const {
        return Versor<T>(this->versor * other.versor);
    }
    
    
    template<typename T>
    c3ga::Mvec<T> Versor<T>::operator()(const c3ga::Mvec<T> &object) const {
        return versor * object * inverse;
    }
    
    
    template<typename T>
    Versor<T> Versor<T>::operator()(const Versor<T> &other) const {
        return Versor<T>((*this)(other.versor));
    }
    
    /////////////////////////// DEDUCTION GUIDES ///////////////////////////////
    
    template<class T>
    Versor(const c3ga::Mvec<T> &) -> Versor<T>;
    
    template<class T>
    Versor(const c3ga::Mvec<T> &, const c3ga::Mvec<T> &) -> Versor<T>;
    
    /////////////////////////////////// ALIAS //////////////////////////////////
    
    using fVersor = Versor<GLfloat>;
    
    using i8Versor = Versor<GLbyte>;
    
    using iu8Versor = Versor<GLubyte>;
    
    using i16Versor = Versor<GLshort>;
    
    using iu16Versor = Versor<GLushort>;
    
    using i32Versor = Versor<GLint>;
    
    using iu32Versor = Versor<GLuint>;
    
    using i64Versor = Versor<GLint64>;
    
    using iu64Versor = Versor<GLuint64>;
}

#endif //GA_PHYSICS_VERSOR_HPP

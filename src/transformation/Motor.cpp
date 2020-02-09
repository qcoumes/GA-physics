#include <transformation/Motor.hpp>


namespace transformation {
    
    template<typename T>
    Motor<T>::Motor(const std::initializer_list<T> &versors) {
        c3ga::Mvec<T> motor = 1;
        
        for (const c3ga::Mvec<T> &versor: versors) {
            motor *= versor;
        }
        
        this->versor = motor;
        this->inverse = motor.inv();
    }
}

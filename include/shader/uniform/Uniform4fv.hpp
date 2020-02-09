#ifndef GA_PHYSICS_UNIFORM4FV_HPP
#define GA_PHYSICS_UNIFORM4FV_HPP

#include <shader/uniform/IUniform.hpp>


namespace mastercraft::shader {
    
    class Uniform4fv : public IUniform {
        
        public:
        
            using IUniform::IUniform;
            
            void load(const void *value) final;
    };
}

#endif //GA_PHYSICS_UNIFORM4FV_HPP

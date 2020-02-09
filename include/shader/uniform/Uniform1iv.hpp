#ifndef GA_PHYSICS_UNIFORM1iv_HPP
#define GA_PHYSICS_UNIFORM1iv_HPP

#include <shader/uniform/IUniform.hpp>


namespace mastercraft::shader {
    
    class Uniform1iv : public IUniform {
        
        public:
        
            using IUniform::IUniform;
            
            void load(const void *value) final;
    };
}

#endif //GA_PHYSICS_UNIFORM1iv_HPP

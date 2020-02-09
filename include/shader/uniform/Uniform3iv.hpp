#ifndef GA_PHYSICS_UNIFORM3iv_HPP
#define GA_PHYSICS_UNIFORM3iv_HPP

#include <shader/uniform/IUniform.hpp>


namespace mastercraft::shader {
    
    class Uniform3iv : public IUniform {
        
        public:
        
            using IUniform::IUniform;
            
            void load(const void *value) final;
    };
}

#endif //GA_PHYSICS_UNIFORM3iv_HPP

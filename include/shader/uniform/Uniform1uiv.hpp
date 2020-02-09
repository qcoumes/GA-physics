#ifndef GA_PHYSICS_UNIFORM1uiv_HPP
#define GA_PHYSICS_UNIFORM1uiv_HPP

#include <shader/uniform/IUniform.hpp>


namespace mastercraft::shader {
    
    class Uniform1uiv : public IUniform {
        
        public:
            
            using IUniform::IUniform;
            
            void load(const void *value) final;
    };
}

#endif //GA_PHYSICS_UNIFORM1uiv_HPP

#ifndef GA_PHYSICS_UNIFORM3uiv_HPP
#define GA_PHYSICS_UNIFORM3uiv_HPP

#include <shader/uniform/IUniform.hpp>


namespace mastercraft::shader {
    
    class Uniform3uiv : public IUniform {
        
        public:
        
            using IUniform::IUniform;
            
            void load(const void *value) final;
    };
}

#endif //GA_PHYSICS_UNIFORM3uiv_HPP

#ifndef GA_PHYSICS_UNIFORMATRIX23_HPP
#define GA_PHYSICS_UNIFORMATRIX23_HPP

#include <shader/uniform/IUniform.hpp>


namespace mastercraft::shader {
    
    class UniformMatrix2x3fv : public IUniform {
        
        public:
        
            using IUniform::IUniform;
            
            void load(const void *value) final;
    };
}

#endif //GA_PHYSICS_UNIFORMATRIX23_HPP

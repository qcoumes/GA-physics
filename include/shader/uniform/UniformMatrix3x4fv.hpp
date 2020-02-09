#ifndef GA_PHYSICS_UNIFORMMATRIX34_HPP
#define GA_PHYSICS_UNIFORMMATRIX34_HPP

#include <shader/uniform/IUniform.hpp>


namespace mastercraft::shader {
    
    class UniformMatrix3x4fv : public IUniform {
        
        public:
        
            using IUniform::IUniform;
            
            void load(const void *value) final;
    };
}

#endif //GA_PHYSICS_UNIFORMMATRIX34_HPP

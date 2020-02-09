#ifndef OPENGL_RENDERED_HPP
#define OPENGL_RENDERED_HPP

#include <app/Versor.hpp>

namespace tool {
    
    class Rendered {
        
        public:
            
            virtual ~Rendered() = default;
            
            virtual void update() = 0;
            
            virtual void render() const = 0;
        
            virtual void transform(const app::Versor<GLfloat> &versor) = 0;
    };
}

#endif //OPENGL_RENDERED_HPP

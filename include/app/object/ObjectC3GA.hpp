#ifndef OPENGL_OBJECTC3GA_HPP
#define OPENGL_OBJECTC3GA_HPP

#include <c3ga/Mvec.hpp>
#include <GL/glew.h>

#include <tool/Rendered.hpp>
#include <app/Versor.hpp>


namespace app::object {
    
    class ObjectC3GA : public tool::Rendered {
        
        public:
        
            virtual void transform(const Versor<GLfloat> &versor) = 0;
        
            [[nodiscard]] virtual c3ga::Mvec<GLfloat> collide(c3ga::Mvec<GLfloat> mvec) = 0;
    };
}

#endif //OPENGL_OBJECTC3GA_HPP

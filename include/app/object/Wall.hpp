#ifndef OPENGL_WALL_HPP
#define OPENGL_WALL_HPP

#include <GL/glew.h>
#include <c3ga/Mvec.hpp>

#include <tool/Rendered.hpp>
#include <shader/ShaderTexture.hpp>
#include <shader/Texture.hpp>


namespace app::object {
    
    class Wall : public tool::Rendered {
        
        private:
            static constexpr GLuint VERTEX_ATTR_POSITION = 0;
            static constexpr GLuint VERTEX_ATTR_NORMAL = 1;
            static constexpr GLuint VERTEX_ATTR_TEXTURE = 2;
            
            shader::ShaderTexture shader;
            shader::Texture texture;
            
            c3ga::Mvec<GLfloat> plane;
            c3ga::Mvec<GLfloat> a;
            c3ga::Mvec<GLfloat> b;
            c3ga::Mvec<GLfloat> c;
            c3ga::Mvec<GLfloat> d;
            
            GLboolean modified;
            GLuint vbo = 0;
            GLuint vao = 0;
        
        public:
            
            Wall();
            
            ~Wall() override;
            
            void transform(const Versor<GLfloat> &versor) override;
            
            void update() override;
            
            void render() const override;
    };
}

#endif //OPENGL_WALL_HPP

#ifndef OPENGL_SQUARE_HPP
#define OPENGL_SQUARE_HPP

#include <GL/glew.h>
#include <c3ga/Mvec.hpp>

#include <shader/ShaderTexture.hpp>
#include <shader/Texture.hpp>
#include <app/object/ObjectC3GA.hpp>


namespace app::object {
    
class Square : public ObjectC3GA {
    
    private:
        static constexpr GLuint VERTEX_ATTR_POSITION = 0;
        static constexpr GLuint VERTEX_ATTR_NORMAL = 1;
        static constexpr GLuint VERTEX_ATTR_TEXTURE = 2;
        
        std::shared_ptr<shader::ShaderTexture> shader;
        std::shared_ptr<shader::Texture> texture;
        
        c3ga::Mvec<GLfloat> plane;
        c3ga::Mvec<GLfloat> normal;
        c3ga::Mvec<GLfloat> a;
        c3ga::Mvec<GLfloat> b;
        c3ga::Mvec<GLfloat> c;
        c3ga::Mvec<GLfloat> d;
        
        GLboolean modified;
        GLuint vbo = 0;
        GLuint vao = 0;
        
        public:
            
            Square(c3ga::Mvec<GLfloat> a, c3ga::Mvec<GLfloat> b, c3ga::Mvec<GLfloat> c, c3ga::Mvec<GLfloat> d,
                   std::shared_ptr<shader::ShaderTexture> shader, std::shared_ptr<shader::Texture> texture);
            
            ~Square() override;
            
            void transform(const Versor<GLfloat> &versor) override;
        
            [[nodiscard]] c3ga::Mvec<GLfloat> collide(c3ga::Mvec<GLfloat> mvec) override;
        
            bool update() override;
            
            void render() const override;
    };
}

#endif //OPENGL_SQUARE_HPP

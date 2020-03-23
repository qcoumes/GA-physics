#ifndef OPENGL_Triangle_HPP
#define OPENGL_Triangle_HPP

#include <GL/glew.h>
#include <c3ga/Mvec.hpp>

#include <shader/ShaderTexture.hpp>
#include <shader/Texture.hpp>
#include <app/object/ObjectC3GA.hpp>


namespace app::object {
    
    class Triangle : public ObjectC3GA {
        
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
            static constexpr GLfloat SIN_60 = 0.8660254037844386467637231707529361834714026269051903140279034897259665f;
        
        public:
            
            Triangle(const c3ga::Mvec<GLfloat> &a, const c3ga::Mvec<GLfloat> &b, const c3ga::Mvec<GLfloat> &c,
                     std::shared_ptr<shader::ShaderTexture> shader, std::shared_ptr<shader::Texture> texture);
            
            ~Triangle() override;
            
            void transform(const Versor<GLfloat> &versor) override;
        
            [[nodiscard]] c3ga::Mvec<GLfloat> collide(c3ga::Mvec<GLfloat> mvec) override;
        
            bool update() override;
            
            void render() const override;
    };
}

#endif //OPENGL_Triangle_HPP

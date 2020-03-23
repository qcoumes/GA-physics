#ifndef GA_PHYSICS_SPHERE_HPP
#define GA_PHYSICS_SPHERE_HPP

#include <vector>

#include <GL/glew.h>
#include <c3ga/Mvec.hpp>

#include <misc/INonCopyable.hpp>
#include <shader/ShaderTexture.hpp>
#include <app/object/ObjectVertex.hpp>
#include <app/object/ObjectC3GA.hpp>


namespace app::object {
    
    class Projectile : public ObjectC3GA {
        
        private:
            static constexpr GLuint VERTEX_ATTR_POSITION = 0;
            static constexpr GLuint VERTEX_ATTR_NORMAL = 1;
            static constexpr GLuint VERTEX_ATTR_TEXTURE = 2;
            
            std::shared_ptr<shader::ShaderTexture> shader;
            std::shared_ptr<shader::Texture> texture;
            
            c3ga::Mvec<GLfloat> multivector;
            GLuint longitude;
            GLuint latitude;
            GLuint64 size;
            GLuint vbo = 0;
            GLuint vao = 0;
            
            fVersor translator;
            GLuint bounce;
        
        private:
        
            bool checkCollision();
            
            void buildVertices(std::vector<ObjectVertex> &vertices) const;
            
            void fillGlData(const std::vector<ObjectVertex> &vertices) const;
        
        public:
            
            Projectile(std::shared_ptr<shader::ShaderTexture> shader, std::shared_ptr<shader::Texture> texture);
            
            ~Projectile() override;
            
            void transform(const Versor<GLfloat> &versor) override;
        
            [[nodiscard]] c3ga::Mvec<GLfloat> collide(c3ga::Mvec<GLfloat> mvec) override;
            
            bool update() override;
            
            void render() const override;
    };
}

#endif // GA_PHYSICS_SPHERE_HPP

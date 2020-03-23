#ifndef OPENGL_ARENA_HPP
#define OPENGL_ARENA_HPP

#include <memory>

#include <GL/glew.h>

#include <shader/ShaderTexture.hpp>
#include <shader/Texture.hpp>
#include <app/object/ObjectC3GA.hpp>
#include <app/Versor.hpp>


namespace app::object {
    
    class Arena : public ObjectC3GA {
        
        private:
            static constexpr GLfloat SQRT_1_2 = 0.70710678118654752440084436210484903928483593768847f;
            static constexpr GLfloat SQRT_1_2d2 = SQRT_1_2 / 2;
            
            std::unordered_map<std::string, std::shared_ptr<ObjectC3GA>> rendered;
            std::shared_ptr<shader::ShaderTexture> shader;
            GLboolean modified;
        
        private:
            
            void xRing(const std::shared_ptr<shader::ShaderTexture> &shader,
                       const std::shared_ptr<shader::Texture> &texture);
            
            void yRing(const std::shared_ptr<shader::ShaderTexture> &shader,
                       const std::shared_ptr<shader::Texture> &texture);
            
            void zRing(const std::shared_ptr<shader::ShaderTexture> &shader,
                       const std::shared_ptr<shader::Texture> &texture);
            
            void fillRings(const std::shared_ptr<shader::ShaderTexture> &shader,
                           const std::shared_ptr<shader::Texture> &texture);
            
            void triangles(const std::shared_ptr<shader::ShaderTexture> &shader,
                           const std::shared_ptr<shader::Texture> &black,
                           const std::shared_ptr<shader::Texture> &white);
        
        public:
            
            Arena(std::shared_ptr<shader::ShaderTexture> shader);
            
            ~Arena() override = default;
            
            void transform(const Versor <GLfloat> &versor) override;
        
            [[nodiscard]] c3ga::Mvec<GLfloat> collide(c3ga::Mvec<GLfloat> mvec) override;
        
            bool update() override;
            
            void render() const override;
    };
}

#endif //OPENGL_ARENA_HPP

#ifndef GA_PHYSICS_SHADERCUBEMAP_HPP
#define GA_PHYSICS_SHADERCUBEMAP_HPP

#include <shader/Shader.hpp>
#include <shader/Cubemap.hpp>


namespace mastercraft::shader {
    
    class ShaderCubemap : public Shader {
        private:
            GLint uTexture = -1;
        
        public:
            
            ShaderCubemap() = default;
            
            ShaderCubemap(const std::string &vsFile, const std::string &fsFile);
            
            void bindCubemap(const Cubemap &texture) const;
            
            void unbindCubemap() const;
    };
}

#endif // GA_PHYSICS_SHADERCUBEMAP_HPP

